/** @file Tracking4D.cpp */
#include "Tracking4D.hpp"
#include "core/ModuleFactory.hpp"
#include "core/detector/GeometryManager.hpp"
#include "core/utils/PlotManager.hpp"
#include "core/utils/Logger.hpp"
#include <future>
#include <map>
#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <sstream>
#include <limits>
#include <chrono>

namespace framework {

    // Cumulative diagnostic counters, only surfaced at DEBUG level in finalize().
    static std::atomic<size_t> stat_seeds_built{0};
    static std::atomic<size_t> stat_seeds_fit_failed{0};
    static std::atomic<size_t> stat_mid_clusters_tested{0};
    static std::atomic<size_t> stat_mid_clusters_passed{0};
    static std::atomic<size_t> stat_tracks_passed_minhits{0};
    static std::atomic<size_t> stat_tracks_failed_finalfit{0};
    static std::atomic<size_t> stat_tracks_killed_by_dedup{0};
    static std::atomic<size_t> stat_tracks_accepted{0};

    // Cumulative wall-clock time per phase of run(), across every batch -
    // module-level timing alone doesn't say WHERE inside Tracking4D the
    // time goes. Only the "parallel fit" phase (the thread_pool->submit()
    // loop below) is genuinely parallelized; bucketing (building the
    // per-event cluster map) and aggregation (collecting futures + filling
    // every plot) both run serially on the main thread, so more threads
    // can't help them. Plain double, not atomic like the counters above:
    // these are only ever touched from run()'s own calling thread, which
    // is never invoked concurrently with itself (Tracking4D doesn't opt
    // into Module::runsConcurrently()).
    static double stat_time_bucketing_s = 0.0;
    static double stat_time_parallel_fit_s = 0.0;
    static double stat_time_aggregate_s = 0.0;
    static double stat_time_save_s = 0.0;

    struct TrackingResult {
        std::vector<std::shared_ptr<Track>> tracks;
        std::vector<int> tracks_per_event;
    };

    Tracking4D::Tracking4D(Configuration cfg, Configuration g_cfg, ThreadPool* pool)
        : Module(std::move(cfg), std::move(g_cfg), pool) {
        
        auto cuts = config.getArray<double>("spatial_cut_abs");
        if (cuts.size() >= 2) {
            spatial_cut_x_ = cuts[0];
            spatial_cut_y_ = cuts[1];
        } else {
            spatial_cut_x_ = config.get<double>("spatial_cut_abs_x", 0.1); 
            spatial_cut_y_ = config.get<double>("spatial_cut_abs_y", 0.1); 
        }
        
        time_cut_abs_ = config.get<double>("time_cut_abs", 690000.0);
        min_hits_ = config.get<size_t>("min_hits_on_track", 5);
        // Accepts either detector names or detector types, as a single
        // bare value or a list (see Configuration::getArray()) - e.g.
        // exclude_detectors = "RD53B_115" or ["caendt5742", "cmsit"].
        exclude_detectors_ = config.getArray<std::string>("exclude_detectors");
        unique_cluster_usage_ = config.get<bool>("unique_cluster_usage", false);

        // Matches corry's own isolated_planes/isolation_cut config: two
        // parallel arrays, zipped by index into a per-detector cut.
        auto isolated_planes = config.getArray<std::string>("isolated_planes");
        auto isolation_cuts = config.getArray<double>("isolation_cut");
        for (size_t i = 0; i < isolated_planes.size() && i < isolation_cuts.size(); ++i) {
            isolation_cut_[isolated_planes[i]] = isolation_cuts[i];
        }

        storage_file_ = config.get<std::string>("storage_file", "");
        save_enabled_ = !storage_file_.empty();
    }

    bool Tracking4D::isExcluded(const DetectorGeo& det) const {
        return std::find(exclude_detectors_.begin(), exclude_detectors_.end(), det.name) != exclude_detectors_.end() ||
               std::find(exclude_detectors_.begin(), exclude_detectors_.end(), det.type) != exclude_detectors_.end();
    }

    Tracking4D::ClusterWriter& Tracking4D::getClusterWriter(const std::string& det_name) {
        auto it = cluster_writers_.find(det_name);
        if (it != cluster_writers_.end()) return it->second;

        H5::Group grp = save_h5_->createGroup("/clusters/" + det_name);
        ClusterWriter w;
        w.event_id  = std::make_unique<AppendableDataset>(grp, "event_id", H5::PredType::NATIVE_UINT64);
        w.cluster_id = std::make_unique<AppendableDataset>(grp, "cluster_id", H5::PredType::NATIVE_UINT64);
        w.col       = std::make_unique<AppendableDataset>(grp, "col", H5::PredType::NATIVE_DOUBLE);
        w.row       = std::make_unique<AppendableDataset>(grp, "row", H5::PredType::NATIVE_DOUBLE);
        w.global_x  = std::make_unique<AppendableDataset>(grp, "global_x", H5::PredType::NATIVE_DOUBLE);
        w.global_y  = std::make_unique<AppendableDataset>(grp, "global_y", H5::PredType::NATIVE_DOUBLE);
        w.global_z  = std::make_unique<AppendableDataset>(grp, "global_z", H5::PredType::NATIVE_DOUBLE);
        w.charge    = std::make_unique<AppendableDataset>(grp, "charge", H5::PredType::NATIVE_DOUBLE);
        w.timestamp = std::make_unique<AppendableDataset>(grp, "timestamp", H5::PredType::NATIVE_DOUBLE);
        // Cluster::seedChannel() - raw CAEN hardware channel of the
        // cluster's highest-charge pixel, needed downstream to determine
        // trigger group membership (see Cluster.hpp's triggerGroupOf()) once
        // reloaded via the Reader module. -1 for non-CAEN detectors/pixels
        // without real channel info.
        w.channel   = std::make_unique<AppendableDataset>(grp, "channel", H5::PredType::NATIVE_INT32);
        auto res = cluster_writers_.emplace(det_name, std::move(w));
        return res.first->second;
    }

    void Tracking4D::saveTracks(const std::vector<std::shared_ptr<Track>>& tracks) {
        if (!save_enabled_ || tracks.empty()) return;

        std::vector<uint64_t> b_event_id, b_track_id;
        std::vector<double> b_mx, b_cx, b_my, b_cy, b_chi2;
        std::vector<int32_t> b_ndof;
        std::vector<uint64_t> b_link_track, b_link_cluster;

        // Per-detector cluster columns, accumulated across every track
        // before a single append() flush per column below - a per-cluster
        // append() would mean one HDF5 extend()+write() call per cluster
        // per column, paying HDF5's per-call metadata overhead on every
        // one; batching removes what would otherwise dominate this
        // module's save-phase cost.
        struct ClusterBuffer {
            std::vector<uint64_t> event_id, cluster_id;
            std::vector<double> col, row, global_x, global_y, global_z, charge, timestamp;
            std::vector<int32_t> channel;
        };
        std::map<std::string, ClusterBuffer> cluster_buffers;

        for (const auto& track : tracks) {
            uint64_t track_id = next_track_id_++;
            b_event_id.push_back(track->eventID());
            b_track_id.push_back(track_id);
            const auto& p = track->params();
            b_mx.push_back(p(0)); b_cx.push_back(p(1)); b_my.push_back(p(2)); b_cy.push_back(p(3));
            b_chi2.push_back(track->chi2());
            b_ndof.push_back(track->ndof());

            for (const auto& cl : track->getClusters()) {
                auto& buf = cluster_buffers[cl->detectorID()];
                uint64_t cluster_id = next_cluster_id_++;
                buf.event_id.push_back(cl->eventID());
                buf.cluster_id.push_back(cluster_id);
                buf.col.push_back(cl->column()); buf.row.push_back(cl->row());
                buf.global_x.push_back(cl->global().x()); buf.global_y.push_back(cl->global().y());
                buf.global_z.push_back(cl->global().z());
                buf.charge.push_back(cl->charge()); buf.timestamp.push_back(cl->timestamp());
                buf.channel.push_back(cl->seedChannel());

                b_link_track.push_back(track_id);
                b_link_cluster.push_back(cluster_id);
            }
        }

        for (auto& [det_name, buf] : cluster_buffers) {
            auto& cw = getClusterWriter(det_name);
            cw.event_id->append(buf.event_id);
            cw.cluster_id->append(buf.cluster_id);
            cw.col->append(buf.col);
            cw.row->append(buf.row);
            cw.global_x->append(buf.global_x);
            cw.global_y->append(buf.global_y);
            cw.global_z->append(buf.global_z);
            cw.charge->append(buf.charge);
            cw.timestamp->append(buf.timestamp);
            cw.channel->append(buf.channel);
        }

        tr_event_id_->append(b_event_id);
        tr_track_id_->append(b_track_id);
        tr_mx_->append(b_mx); tr_cx_->append(b_cx); tr_my_->append(b_my); tr_cy_->append(b_cy);
        tr_chi2_->append(b_chi2);
        tr_ndof_->append(b_ndof);
        link_track_id_->append(b_link_track);
        link_cluster_id_->append(b_link_cluster);
    }

    void Tracking4D::saveDutClusters(const DataBatch& batch) {
        if (!save_enabled_) return;
        auto& geo = GeometryManager::getInstance();

        for (auto const& [det_name, clusters] : batch.clusters) {
            if (clusters.empty() || !geo.hasDetector(det_name)) continue;
            auto& det = geo.getDetector(det_name);
            // Any detector excluded from the fit (via #exclude_detectors_)
            // whose clusters ended up in batch.clusters anyway - includes
            // CAEN boards now that WaveformSelector feeds their selected
            // waveforms into ClusteringSpatial (as Pixels, upstream of this
            // module) instead of building clusters directly.
            if (!isExcluded(det)) continue;

            std::vector<uint64_t> b_ev, b_cid;
            std::vector<double> b_col, b_row, b_gx, b_gy, b_gz, b_chg, b_ts;
            std::vector<int32_t> b_channel;
            for (auto const& cl : clusters) {
                b_ev.push_back(cl->eventID());
                b_cid.push_back(next_cluster_id_++);
                b_col.push_back(cl->column()); b_row.push_back(cl->row());
                b_gx.push_back(cl->global().x()); b_gy.push_back(cl->global().y()); b_gz.push_back(cl->global().z());
                b_chg.push_back(cl->charge()); b_ts.push_back(cl->timestamp());
                b_channel.push_back(cl->seedChannel());
            }

            auto& cw = getClusterWriter(det_name);
            cw.event_id->append(b_ev); cw.cluster_id->append(b_cid);
            cw.col->append(b_col); cw.row->append(b_row);
            cw.global_x->append(b_gx); cw.global_y->append(b_gy); cw.global_z->append(b_gz);
            cw.charge->append(b_chg); cw.timestamp->append(b_ts);
            cw.channel->append(b_channel);
        }
    }

    void Tracking4D::initialize() {
        WR_LOG(STATUS, "Initializing Tracking4D...");
        auto& pm = PlotManager::getInstance();
        pm.registerPlot1D(getName(), "track_chi2", 300, 0, 150);
        pm.registerPlot1D(getName(), "track_chi2ndof", 500, 0, 50);
        pm.registerPlot1D(getName(), "clusters_per_track", 10, -0.5, 9.5);
        pm.registerPlot1D(getName(), "tracks_per_event", 100, -0.5, 99.5);
        pm.registerPlot1D(getName(), "track_angle_x", 2000, -0.01, 0.01);
        pm.registerPlot1D(getName(), "track_angle_y", 2000, -0.01, 0.01);
        // Matches corry's tracksVsTime - unlike trackTime/trackTimeTrigger/
        // trackTimeTriggerChi2 (also in corry's Tracking4D), this doesn't
        // need an Event/trigger-list concept Warlock doesn't have yet -
        // just the track's own timestamp, same as every other timestamp
        // plot in this codebase.
        pm.registerPlot1D(getName(), "tracks_vs_time", 3000000, 0, 3000);

        auto& geo = GeometryManager::getInstance();
        for (const auto& name : geo.getDetectorNames()) {
            auto& det = geo.getDetector(name);
            if (!det.isPixelDetector()) continue;
            std::string dir = getName() + "/" + name;

            // Matches corry's own scope: local_intersect/global_intersect/
            // isolation_distance are registered for every regular pixel
            // detector, including ones excluded from the fit (DUTs) -
            // corry's get_regular_detectors(true) - so you can see where a
            // fitted track (built from other planes) would land on a DUT
            // that never contributed to the fit. Only the residual/pull
            // family below is gated to non-excluded (fit-participating)
            // planes.
            pm.registerPlot2D(dir, "local_intersect", det.n_pixels_x, 0, det.n_pixels_x, det.n_pixels_y, 0, det.n_pixels_y);
            pm.registerPlot2D(dir, "global_intersect", 600, -30, 30, 600, -30, 30);
            // corry's IsolationDistance is registered unconditionally too
            // (only ever FILLED if isolation_cut_ is non-empty - see run()).
            pm.registerPlot1D(dir, "isolation_distance", 500, 0, 5000);

            if (isExcluded(det)) continue;

            // Matches corry's own GlobalResidualsX/Y (and Local, same
            // range) exactly (Tracking4D.cpp: `3 * detector->getPitch().X()/Y()`).
            double range_x = det.pitch_x * 3.0;
            double range_y = det.pitch_y * 3.0;

            pm.registerPlot1D(dir, "residual_x_global", 500, -range_x, range_x);
            pm.registerPlot1D(dir, "residual_y_global", 500, -range_y, range_y);
            pm.registerPlot1D(dir, "residual_z_global", 500, -0.1, 0.1);
            pm.registerPlot2D(dir, "residual_x_vs_position_x_global", 500, -range_x, range_x, 400, -10.0, 10.0);
            pm.registerPlot2D(dir, "residual_y_vs_position_y_global", 500, -range_y, range_y, 400, -10.0, 10.0);
            // Cross-axis: X residual vs Y position, Y residual vs X position -
            // corry has both same-axis AND cross-axis views.
            pm.registerPlot2D(dir, "residual_x_vs_position_y_global", 500, -range_x, range_x, 400, -10.0, 10.0);
            pm.registerPlot2D(dir, "residual_y_vs_position_x_global", 500, -range_y, range_y, 400, -10.0, 10.0);

            // Local residuals - the detector-frame equivalent (rotation
            // removed) of the global ones above; the two only coincide
            // exactly for an unrotated detector.
            pm.registerPlot1D(dir, "residual_x_local", 500, -range_x, range_x);
            pm.registerPlot1D(dir, "residual_y_local", 500, -range_y, range_y);

            // Cluster-width breakdown, both frames - matches corry's
            // residualsXwidth1/2/3 (and Y), local and global.
            for (int w = 1; w <= 3; ++w) {
                pm.registerPlot1D(dir, "residual_x_local_width" + std::to_string(w), 500, -range_x, range_x);
                pm.registerPlot1D(dir, "residual_y_local_width" + std::to_string(w), 500, -range_y, range_y);
                pm.registerPlot1D(dir, "residual_x_global_width" + std::to_string(w), 500, -range_x, range_x);
                pm.registerPlot1D(dir, "residual_y_global_width" + std::to_string(w), 500, -range_y, range_y);
            }

            // Pulls - residual divided by the cluster's own position error.
            // Warlock has no per-cluster error (corry's can vary with
            // column/row for some detector types); DetectorGeo::spatial_res_x/y
            // (the detector's configured nominal resolution) stands in for
            // it here, which is what corry's own getSpatialResolution()
            // reduces to for a fixed-resolution detector anyway.
            pm.registerPlot1D(dir, "pull_x_local", 500, -5, 5);
            pm.registerPlot1D(dir, "pull_y_local", 500, -5, 5);
            pm.registerPlot1D(dir, "pull_x_global", 500, -5, 5);
            pm.registerPlot1D(dir, "pull_y_global", 500, -5, 5);
        }

        stat_seeds_built = 0; stat_seeds_fit_failed = 0; stat_mid_clusters_tested = 0;
        stat_mid_clusters_passed = 0; stat_tracks_passed_minhits = 0; stat_tracks_failed_finalfit = 0;
        stat_tracks_killed_by_dedup = 0; stat_tracks_accepted = 0;

        if (save_enabled_) {
            save_h5_ = std::make_unique<H5::H5File>(storage_file_, H5F_ACC_TRUNC);
            H5::Group tr_grp = save_h5_->createGroup("/tracks");
            tr_event_id_ = std::make_unique<AppendableDataset>(tr_grp, "event_id", H5::PredType::NATIVE_UINT64);
            tr_track_id_ = std::make_unique<AppendableDataset>(tr_grp, "track_id", H5::PredType::NATIVE_UINT64);
            tr_mx_   = std::make_unique<AppendableDataset>(tr_grp, "mx", H5::PredType::NATIVE_DOUBLE);
            tr_cx_   = std::make_unique<AppendableDataset>(tr_grp, "cx", H5::PredType::NATIVE_DOUBLE);
            tr_my_   = std::make_unique<AppendableDataset>(tr_grp, "my", H5::PredType::NATIVE_DOUBLE);
            tr_cy_   = std::make_unique<AppendableDataset>(tr_grp, "cy", H5::PredType::NATIVE_DOUBLE);
            tr_chi2_ = std::make_unique<AppendableDataset>(tr_grp, "chi2", H5::PredType::NATIVE_DOUBLE);
            tr_ndof_ = std::make_unique<AppendableDataset>(tr_grp, "ndof", H5::PredType::NATIVE_INT32);

            H5::Group link_grp = save_h5_->createGroup("/track_cluster_links");
            link_track_id_   = std::make_unique<AppendableDataset>(link_grp, "track_id", H5::PredType::NATIVE_UINT64);
            link_cluster_id_ = std::make_unique<AppendableDataset>(link_grp, "cluster_id", H5::PredType::NATIVE_UINT64);

            // Must exist before getClusterWriter() can create "/clusters/<det_name>"
            // groups under it - createGroup() doesn't auto-create missing
            // parent groups.
            save_h5_->createGroup("/clusters");

            WR_LOG(STATUS, "Tracking4D: saving fitted tracks to " + storage_file_);
        }
    }

    void Tracking4D::run(DataBatch& batch) {
        auto t_bucket_start = std::chrono::steady_clock::now();

        auto& geo = GeometryManager::getInstance();
        std::vector<std::string> planes;
        // Every regular pixel detector, INCLUDING excluded/DUT ones - used
        // only for post-fit intercept evaluation (local_intersect/
        // global_intersect/isolation_distance), never for building the fit
        // itself. Matches corry's get_regular_detectors(true) scope for
        // those specific plots - see initialize()'s docs.
        std::vector<std::string> all_planes;
        for (auto const& n : geo.getDetectorNames()) {
            auto& det = geo.getDetector(n);
            if (!det.isPixelDetector()) continue;
            all_planes.push_back(n);
            if (!isExcluded(det)) planes.push_back(n);
        }
        if (planes.size() < 2) return;

        std::sort(planes.begin(), planes.end(), [&](const std::string& a, const std::string& b) {
            return geo.getDetector(a).position.z() < geo.getDetector(b).position.z();
        });
        std::sort(all_planes.begin(), all_planes.end(), [&](const std::string& a, const std::string& b) {
            return geo.getDetector(a).position.z() < geo.getDetector(b).position.z();
        });

        using EventMap = std::unordered_map<std::string, std::vector<std::shared_ptr<Cluster>>>;
        std::unordered_map<uint64_t, EventMap> events;
        for (auto const& entry : batch.clusters) {
            for (auto const& cl : entry.second) events[cl->eventID()][entry.first].push_back(cl);
        }

        std::vector<std::pair<uint64_t, const EventMap*>> event_list;
        event_list.reserve(events.size());
        for (auto const& entry : events) event_list.push_back({entry.first, &entry.second});

        stat_time_bucketing_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_bucket_start).count();

        size_t n_threads = thread_pool->getThreadCount();
        size_t chunk_size = (event_list.size() + n_threads - 1) / n_threads;

        std::vector<std::future<TrackingResult>> futures;

        for (size_t i = 0; i < event_list.size(); i += chunk_size) {
            size_t end = std::min(i + chunk_size, event_list.size());
            
            futures.push_back(thread_pool->submit([this, i, end, &event_list, &planes, &all_planes]() {
                TrackingResult local_res;
                
                size_t l_seeds_built = 0;
                size_t l_seeds_fit_failed = 0;
                size_t l_mid_tested = 0;
                size_t l_mid_passed = 0;
                size_t l_passed_minhits = 0;
                size_t l_failed_finalfit = 0;
                size_t l_killed_by_dedup = 0;
                size_t l_accepted = 0;

                for (size_t ev_idx = i; ev_idx < end; ++ev_idx) {
                    uint64_t ev_id = event_list[ev_idx].first;
                    const auto& dets = *(event_list[ev_idx].second);

                    std::vector<std::shared_ptr<Track>> candidate_tracks;

                    // Determine the reference pair: the first and last plane (in z-order)
                    // that actually carry clusters in this event. Corryvreckan seeds tracks
                    // only from this single outermost pair, not every possible pair of
                    // planes, to match its own track-finding behaviour.
                    size_t f_idx = std::numeric_limits<size_t>::max(), l_idx = std::numeric_limits<size_t>::max();
                    for (size_t idx = 0; idx < planes.size(); ++idx) {
                        if (dets.count(planes[idx]) == 0) continue;
                        if (f_idx == std::numeric_limits<size_t>::max()) f_idx = idx;
                        l_idx = idx;
                    }

                    if (f_idx != std::numeric_limits<size_t>::max() && f_idx != l_idx && (l_idx - f_idx + 1) >= min_hits_) {
                        std::string first = planes[f_idx];
                        std::string last = planes[l_idx];

                        for (auto& c_first : dets.at(first)) {
                            for (auto& c_last : dets.at(last)) {

                                // Timing cut on the reference pair, so clusters from
                                // unrelated times can't be combined into a seed.
                                if (std::abs(c_first->timestamp() - c_last->timestamp()) > time_cut_abs_) continue;

                                l_seeds_built++;
                                Track candidate;
                                candidate.addCluster(c_first);
                                candidate.addCluster(c_last);
                                candidate.setEventID(ev_id);
                                double ts_sum = c_first->timestamp() + c_last->timestamp();
                                size_t ts_count = 2;
                                candidate.setTimestamp(ts_sum / static_cast<double>(ts_count));

                                if (!candidate.fit()) {
                                    l_seeds_fit_failed++;
                                    continue;
                                }

                                for (size_t mid_idx = f_idx + 1; mid_idx < l_idx; ++mid_idx) {
                                    const std::string& p_name = planes[mid_idx];
                                    if (dets.count(p_name) == 0) continue;

                                    auto* det = dets.at(p_name).front()->detector();
                                    Eigen::Vector3d global_pred = candidate.interceptAt(*det);

                                    if (std::isnan(global_pred.x()) || std::isnan(global_pred.y())) continue;

                                    Eigen::Vector3d local_pred = det->R_inv * (global_pred - det->position);

                                    std::shared_ptr<Cluster> best_match = nullptr;
                                    double closestClusterDistanceSq = spatial_cut_x_*spatial_cut_x_ + spatial_cut_y_*spatial_cut_y_;
                                    double best_ts = 0.0;
                                    double running_ts = ts_sum / static_cast<double>(ts_count);

                                    for (auto& c_mid : dets.at(p_name)) {
                                        l_mid_tested++;

                                        if (std::abs(c_mid->timestamp() - running_ts) > time_cut_abs_) {
                                            continue;
                                        }

                                        Eigen::Vector3d local_hit = det->R_inv * (c_mid->global() - det->position);

                                        double dx = local_hit.x() - local_pred.x();
                                        double dy = local_hit.y() - local_pred.y();

                                        // Elliptical spatial cut, matching corryvreckan:
                                        // (dx/cutx)^2 + (dy/cuty)^2 <= 1 is inside the cut -
                                        // stricter at the corners than a plain box cut.
                                        double norm = (dx * dx) / (spatial_cut_x_ * spatial_cut_x_) +
                                                      (dy * dy) / (spatial_cut_y_ * spatial_cut_y_);
                                        if (norm > 1.0) {
                                            continue;
                                        }

                                        l_mid_passed++;

                                        double distSq = dx * dx + dy * dy;
                                        if (distSq < closestClusterDistanceSq) {
                                            closestClusterDistanceSq = distSq;
                                            best_match = c_mid;
                                            best_ts = c_mid->timestamp();
                                        }
                                    }
                                    if (best_match) {
                                        candidate.addCluster(best_match);
                                        ts_sum += best_ts;
                                        ts_count++;
                                        candidate.setTimestamp(ts_sum / static_cast<double>(ts_count));
                                        // corry refits its reference trajectory after every mid-plane
                                        // match (Track::updatePlane() calls fit() internally) - so
                                        // later mid-plane predictions are based on an increasingly
                                        // accurate line, not the original 2-cluster seed. Return value
                                        // is intentionally ignored here, matching corry's own void
                                        // updatePlane(): any fundamental failure is still caught by the
                                        // final candidate.fit() check below.
                                        candidate.fit();
                                    }
                                }

                                if (candidate.getClusters().size() >= min_hits_) {
                                    l_passed_minhits++;

                                    if (candidate.fit()) {
                                        candidate_tracks.push_back(std::make_shared<Track>(candidate));
                                    } else {
                                        l_failed_finalfit++;
                                    }
                                }
                            }
                        }
                    }

                    int accepted_for_this_event = 0;
                    if (!candidate_tracks.empty()) {
                        // Corryvreckan only rejects tracks sharing clusters when
                        // `unique_cluster_usage` is explicitly enabled (default: off) -
                        // running this unconditionally would silently drop legitimate
                        // overlapping tracks corry keeps by default.
                        if (unique_cluster_usage_ && candidate_tracks.size() > 1) {
                            std::sort(candidate_tracks.begin(), candidate_tracks.end(),
                                      [](const std::shared_ptr<Track>& a, const std::shared_ptr<Track>& b) {
                                          if (a->getClusters().size() != b->getClusters().size()) {
                                              return a->getClusters().size() > b->getClusters().size();
                                          }
                                          return a->getChi2ndof() < b->getChi2ndof();
                                      });

                            std::unordered_set<Cluster*> used_clusters;
                            for (auto& track : candidate_tracks) {
                                bool is_shared = false;
                                for (const auto& c : track->getClusters()) {
                                    if (used_clusters.count(c.get())) {
                                        is_shared = true;
                                        break;
                                    }
                                }

                                if (!is_shared) {
                                    for (const auto& c : track->getClusters()) used_clusters.insert(c.get());
                                    local_res.tracks.push_back(std::move(track));
                                    accepted_for_this_event++;
                                    l_accepted++;
                                } else {
                                    l_killed_by_dedup++;
                                }
                            }
                        } else {
                            for (auto& track : candidate_tracks) {
                                local_res.tracks.push_back(std::move(track));
                                accepted_for_this_event++;
                                l_accepted++;
                            }
                        }
                    }

                    // Track isolation: pairwise distance between every pair
                    // of tracks accepted for THIS event (never cross-event -
                    // matches corry, which evaluates this over its own
                    // per-event `tracks` vector before moving to the next
                    // event), at each configured isolated_planes_ detector.
                    // Only runs at all if isolation_cut_ is non-empty,
                    // matching corry's own "isolation_cut_.size() > 0" gate.
                    if (!isolation_cut_.empty() && accepted_for_this_event > 1) {
                        auto& pm_local = PlotManager::getInstance();
                        auto& geo_local = GeometryManager::getInstance();
                        size_t first_idx = local_res.tracks.size() - static_cast<size_t>(accepted_for_this_event);
                        for (size_t i1 = first_idx; i1 < local_res.tracks.size(); ++i1) {
                            for (size_t i2 = i1 + 1; i2 < local_res.tracks.size(); ++i2) {
                                auto& t1 = local_res.tracks[i1];
                                auto& t2 = local_res.tracks[i2];
                                for (auto const& [det_name, cut] : isolation_cut_) {
                                    if (!geo_local.hasDetector(det_name)) continue;
                                    auto& det = geo_local.getDetector(det_name);
                                    Eigen::Vector3d p1 = t1->interceptAt(det);
                                    Eigen::Vector3d p2 = t2->interceptAt(det);
                                    if (std::isnan(p1.x()) || std::isnan(p2.x())) continue;
                                    Eigen::Vector3d l1 = det.R_inv * (p1 - det.position);
                                    Eigen::Vector3d l2 = det.R_inv * (p2 - det.position);
                                    double dx = l1.x() - l2.x(), dy = l1.y() - l2.y();
                                    double distance_um = std::sqrt(dx * dx + dy * dy) * 1000.0;
                                    pm_local.fill1D(getName() + "/" + det_name, "isolation_distance", distance_um);
                                    if (distance_um < cut) {
                                        t1->markAsNotIsolatedAtPlane(det_name);
                                        t2->markAsNotIsolatedAtPlane(det_name);
                                    }
                                }
                            }
                        }
                    }

                    local_res.tracks_per_event.push_back(accepted_for_this_event);
                }

                stat_seeds_built += l_seeds_built;
                stat_seeds_fit_failed += l_seeds_fit_failed;
                stat_mid_clusters_tested += l_mid_tested;
                stat_mid_clusters_passed += l_mid_passed;
                stat_tracks_passed_minhits += l_passed_minhits;
                stat_tracks_failed_finalfit += l_failed_finalfit;
                stat_tracks_killed_by_dedup += l_killed_by_dedup;
                stat_tracks_accepted += l_accepted;

                // Plot fills happen here, inside the worker, rather than in
                // a later serial pass: PlotManager locks per-histogram (see
                // fill1D/fill2D), so concurrent fills from different
                // chunks are already safe, and filling here overlaps this
                // cost with other chunks' fit work instead of adding pure
                // serial time afterward. Tracks failing the intercept
                // check are dropped here too - batch.tracks only ever
                // wanted the survivors.
                auto& pm = PlotManager::getInstance();
                auto& geo_ref = GeometryManager::getInstance();
                std::vector<std::shared_ptr<Track>> kept_tracks;
                kept_tracks.reserve(local_res.tracks.size());

                for (int n : local_res.tracks_per_event) {
                    pm.fill1D(getName(), "tracks_per_event", static_cast<double>(n));
                }

                for (const auto& track : local_res.tracks) {
                    Eigen::Vector3d p0 = track->interceptAt(geo_ref.getDetector(planes[0]));
                    if (std::isnan(p0.x()) || std::isnan(p0.y())) continue;
                    kept_tracks.push_back(track);

                    pm.fill1D(getName(), "track_chi2", track->chi2());
                    pm.fill1D(getName(), "track_chi2ndof", track->getChi2ndof());
                    pm.fill1D(getName(), "clusters_per_track", static_cast<double>(track->getClusters().size()));
                    pm.fill1D(getName(), "track_angle_x", std::atan(track->params()(0)));
                    pm.fill1D(getName(), "track_angle_y", std::atan(track->params()(2)));
                    pm.fill1D(getName(), "tracks_vs_time", track->timestamp() / 1.0e9);

                    // Per-detector residuals/pulls, both frames (matches
                    // corryvreckan's residualsX/Y_local/_global, the
                    // cluster-width breakdowns, and the pulls). Only ever
                    // touches fit-participating (non-excluded) detectors,
                    // since track->getClusters() by construction only ever
                    // contains clusters from `planes`, never `all_planes`'
                    // excluded members - matches the registration gating in
                    // initialize().
                    for (const auto& cl : track->getClusters()) {
                        auto* det = cl->detector();
                        std::string dir = getName() + "/" + det->name;
                        Eigen::Vector3d pred = track->interceptAt(*det);
                        double res_x = cl->global().x() - pred.x();
                        double res_y = cl->global().y() - pred.y();
                        double res_z = cl->global().z() - pred.z();
                        pm.fill1D(dir, "residual_x_global", res_x);
                        pm.fill1D(dir, "residual_y_global", res_y);
                        pm.fill1D(dir, "residual_z_global", res_z);
                        pm.fill2D(dir, "residual_x_vs_position_x_global", res_x, cl->global().x());
                        pm.fill2D(dir, "residual_y_vs_position_y_global", res_y, cl->global().y());
                        pm.fill2D(dir, "residual_x_vs_position_y_global", res_x, cl->global().y());
                        pm.fill2D(dir, "residual_y_vs_position_x_global", res_y, cl->global().x());

                        // Local frame: same residual, with the detector's
                        // own rotation removed - coincides with the global
                        // one exactly only for an unrotated detector.
                        Eigen::Vector3d local_pred = det->R_inv * (pred - det->position);
                        Eigen::Vector3d local_hit = det->R_inv * (cl->global() - det->position);
                        double res_x_local = local_hit.x() - local_pred.x();
                        double res_y_local = local_hit.y() - local_pred.y();
                        pm.fill1D(dir, "residual_x_local", res_x_local);
                        pm.fill1D(dir, "residual_y_local", res_y_local);

                        // Pulls: residual over the detector's configured
                        // spatial resolution - see initialize()'s docs on
                        // why this stands in for corry's per-cluster error.
                        pm.fill1D(dir, "pull_x_local", res_x_local / det->spatial_res_x);
                        pm.fill1D(dir, "pull_y_local", res_y_local / det->spatial_res_y);
                        pm.fill1D(dir, "pull_x_global", res_x / det->spatial_res_x);
                        pm.fill1D(dir, "pull_y_global", res_y / det->spatial_res_y);

                        int col_width = cl->columnWidth();
                        int row_width = cl->rowWidth();
                        if (col_width >= 1 && col_width <= 3) {
                            pm.fill1D(dir, "residual_x_local_width" + std::to_string(col_width), res_x_local);
                            pm.fill1D(dir, "residual_x_global_width" + std::to_string(col_width), res_x);
                        }
                        if (row_width >= 1 && row_width <= 3) {
                            pm.fill1D(dir, "residual_y_local_width" + std::to_string(row_width), res_y_local);
                            pm.fill1D(dir, "residual_y_global_width" + std::to_string(row_width), res_y);
                        }
                    }

                    // Track intercepts on every regular pixel detector
                    // (including excluded/DUT ones - see all_planes' docs),
                    // whether or not the track has a cluster there, local
                    // col/row and global XY.
                    for (const auto& p_name : all_planes) {
                        auto& det = geo_ref.getDetector(p_name);
                        Eigen::Vector3d pred = track->interceptAt(det);
                        if (std::isnan(pred.x()) || std::isnan(pred.y())) continue;

                        std::string dir = getName() + "/" + p_name;
                        Eigen::Vector3d local_pred = det.R_inv * (pred - det.position);
                        double col = det.getColumn(local_pred.x());
                        double row = det.getRow(local_pred.y());
                        pm.fill2D(dir, "local_intersect", col, row);
                        pm.fill2D(dir, "global_intersect", pred.x(), pred.y());
                    }
                }
                local_res.tracks = std::move(kept_tracks);
                total_tracks_ += local_res.tracks.size();

                return local_res;
            }));
        }

        // Waiting for every chunk's fit result is the parallelized phase's
        // real wall-clock cost - kept separate from the plotting loop
        // below so each can be timed on its own. Collecting all results
        // before processing any of them doesn't change what gets plotted,
        // or in what order, since chunks are independent.
        auto t_wait_start = std::chrono::steady_clock::now();
        std::vector<TrackingResult> results;
        results.reserve(futures.size());
        for (auto& f : futures) results.push_back(f.get());
        stat_time_parallel_fit_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_wait_start).count();

        // Plot fills happen inside each worker (see above), so all that's
        // left here is merging every chunk's surviving tracks into the
        // batch - cheap pointer-vector concatenation, not per-track
        // computation.
        auto t_agg_start = std::chrono::steady_clock::now();
        for (auto& result : results) {
            batch.tracks.insert(batch.tracks.end(), result.tracks.begin(), result.tracks.end());
        }
        stat_time_aggregate_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_agg_start).count();

        auto t_save_start = std::chrono::steady_clock::now();
        if (save_enabled_) {
            saveTracks(batch.tracks);
            saveDutClusters(batch);
        }
        stat_time_save_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_save_start).count();
    }

    void Tracking4D::finalize() {
        std::ostringstream stats;
        stats << "Pipeline summary: "
              << "seeds_built=" << stat_seeds_built.load()
              << " seeds_failed_fit=" << stat_seeds_fit_failed.load()
              << " mid_clusters_tested=" << stat_mid_clusters_tested.load()
              << " mid_clusters_passed=" << stat_mid_clusters_passed.load()
              << " tracks_passed_minhits=" << stat_tracks_passed_minhits.load()
              << " tracks_failed_finalfit=" << stat_tracks_failed_finalfit.load()
              << " tracks_killed_by_dedup=" << stat_tracks_killed_by_dedup.load()
              << " tracks_accepted=" << stat_tracks_accepted.load();
        WR_LOG(DEBUG, stats.str());

        WR_LOG(STATUS, "Tracking finished. Total tracks found: " + std::to_string(total_tracks_));

        std::ostringstream phase_stats;
        phase_stats << "Tracking4D phase timing (seconds): "
                    << "bucketing=" << stat_time_bucketing_s
                    << " parallel_fit_and_plot_wait=" << stat_time_parallel_fit_s
                    << " aggregate_merge=" << stat_time_aggregate_s
                    << " save=" << stat_time_save_s;
        WR_LOG(STATUS, phase_stats.str());

        if (save_enabled_ && save_h5_) {
            save_h5_->close();
            WR_LOG(STATUS, "Tracking4D: saved " + std::to_string(next_track_id_) + " tracks (" +
                        std::to_string(next_cluster_id_) + " cluster rows) to " + storage_file_);
        }
    }

    REGISTER_MODULE(Tracking4D)
}