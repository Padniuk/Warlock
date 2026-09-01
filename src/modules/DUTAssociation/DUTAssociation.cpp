#include "DUTAssociation.hpp"
#include "core/ModuleFactory.hpp"
#include "core/detector/GeometryManager.hpp"
#include "core/utils/PlotManager.hpp"
#include "core/utils/Logger.hpp"
#include <future>
#include <limits>
#include <cmath>
#include <set>

namespace framework {
    DUTAssociation::DUTAssociation(Configuration cfg, Configuration g_cfg, ThreadPool* pool)
        : Module(std::move(cfg), std::move(g_cfg), pool) {

        dut_name_ = config.get<std::string>("name", "");
        type_filter_ = config.get<std::string>("type", "");

        auto cuts = config.getArray<double>("spatial_cut_abs");
        spatial_cut_x_ = (cuts.size() > 0) ? cuts[0] : 1.0;
        spatial_cut_y_ = (cuts.size() > 1) ? cuts[1] : 1.0;

        // corry defaults this cut effectively off (1e30) unless configured.
        time_cut_abs_ = config.get<double>("time_cut_abs", 1e30);
        use_cluster_centre_ = config.get<bool>("use_cluster_centre", false);
    }

    void DUTAssociation::initialize() {
        auto& geo = GeometryManager::getInstance();

        // Bound either to exactly one detector ("name") or auto-instantiated
        // once per matching detector ("type") - same resolution as
        // AnalysisEfficiency, matching corry's own module-instantiation
        // semantics (ModuleManager.cpp).
        if (!dut_name_.empty()) {
            if (!geo.hasDetector(dut_name_)) {
                throw std::runtime_error("DUTAssociation: detector '" + dut_name_ + "' not found in geometry.");
            }
            target_detectors_.push_back(dut_name_);
        } else if (!type_filter_.empty()) {
            for (const auto& name : geo.getDetectorNames()) {
                if (geo.getDetector(name).type == type_filter_) target_detectors_.push_back(name);
            }
        } else {
            throw std::runtime_error("DUTAssociation: must configure either 'name' or 'type'.");
        }
        WR_LOG(STATUS, "Initializing Association for: " + (dut_name_.empty() ? "type " + type_filter_ : dut_name_));

        for (const auto& name : target_detectors_) {
            auto& det = geo.getDetector(name);
            auto& pm = PlotManager::getInstance();

            // Registers the full histogram set for one plot subpath -
            // called once for the combined (both-dies) view and, for
            // two-die boards (DetectorGeo::dieOf()), once more each for
            // "<name>/TOP" and "<name>/BOTTOM" - see run() for how the
            // per-cluster fills are routed to these.
            auto registerFor = [&](const std::string& dir) {
                // Binning matched to corry's core DUTAssociation histograms.
                pm.registerPlot1D(dir, "hDistXClusterClosestPx", 2000, -1000, 1000);
                pm.registerPlot1D(dir, "hDistYClusterClosestPx", 2000, -1000, 1000);
                pm.registerPlot1D(dir, "no_assoc_cls", 10, -0.5, 9.5);
                pm.registerPlot1D(dir, "hCutHisto", 2, 1, 3);

                // 2D local track-cluster distance, filled for every examined
                // cluster regardless of cut outcome (corry's "before cuts" view),
                // at two different zoom levels.
                pm.registerPlot2D(dir, "hDist_trackCluster_2D", 500, -det.pitch_x, det.pitch_x, 500, -det.pitch_x, det.pitch_x);
                pm.registerPlot2D(dir, "hDist_trackCluster_2D_huge",
                                   500, -det.pitch_x * 5, det.pitch_x * 5, 500, -det.pitch_x * 5, det.pitch_x * 5);
                pm.registerPlot2D(dir, "hDist_trackCluster_2D_huge_globL",
                                   500, -det.pitch_x * 5, det.pitch_x * 5, 500, -det.pitch_x * 5, det.pitch_x * 5);

                // Filled only for successfully associated (track, cluster) pairs.
                pm.registerPlot2D(dir, "hDist_trackCluster_2D_assoc", 500, -det.pitch_x, det.pitch_x, 500, -det.pitch_x, det.pitch_x);
                pm.registerPlot2D(dir, "hDist_trackCluster_2D_assoc_global", 500, -det.pitch_x, det.pitch_x, 500, -det.pitch_x, det.pitch_x);
                double size_x = det.sizeX(), size_y = det.sizeY();
                pm.registerPlot2D(dir, "hClusterPosAssoc2D", 600, -size_x * 1.1, size_x * 1.1, 600, -size_y * 1.1, size_y * 1.1);
                // Matches corry's own DUTAssociation.cpp exactly: both
                // hTrackPosAssoc2D and hTrackPosNoAssoc2D use the same fixed
                // generic +-10mm range (only hClusterPosAssoc2D is sized to
                // the board, via getSize()*1.1 there too).
                pm.registerPlot2D(dir, "hTrackPosAssoc2D", 600, -10, 10, 600, -10, 10);
                pm.registerPlot2D(dir, "hTrackPosNoAssoc2D", 600, -10, 10, 600, -10, 10);
            };

            std::string dir = getName() + "/" + name;
            registerFor(dir);
            // Single shared source of truth (DetectorGeo::dieLabels(), see
            // GeometryManager.hpp) instead of re-deriving this
            // detector-by-detector - AnalysisEfficiency/AnalysisTiming
            // check the exact same thing.
            auto labels = det.dieLabels();
            // Only warn for a genuinely unrecognized layout string - not
            // for a TILGAD/TREF board that's simply single-row (e.g.
            // CAEN_UZH_2's genuinely single-row TILGAD board), where
            // dieLabels() correctly returns empty with nothing wrong.
            if (det.hasUnrecognizedLayout()) {
                WR_LOG(WARNING, "DUTAssociation: detector '" + name + "' has unrecognized layout '" + det.layout +
                                     "' - not splitting into per-die plots. Add a dedicated Layout subclass "
                                     "(src/core/detector/layouts/) if this layout needs die-level splitting.");
            }
            for (auto const& die : labels) registerFor(dir + "/" + die);
            die_labels_by_det_[name] = std::move(labels);

            WR_LOG(DEBUG, "Monitoring enabled for DUT: " + name);
        }
    }

    void DUTAssociation::run(DataBatch& batch) {
        if (batch.tracks.empty()) return;

        auto& geo = GeometryManager::getInstance();
        std::vector<std::string> target_detectors;
        for (const auto& name : target_detectors_) {
            if (batch.clusters.count(name) > 0) target_detectors.push_back(name);
        }
        if (target_detectors.empty()) return;

        // Group each target detector's clusters by event once per batch - a
        // track must only ever be compared against clusters from the same
        // event. Without this, every track would be compared against every
        // cluster in the whole batch, an O(n_tracks * n_clusters) blowup
        // rejected only afterwards by the timing cut, unlike corry's
        // per-event clipboard model which never hits this by construction.
        std::map<std::string, std::map<uint64_t, std::vector<std::shared_ptr<Cluster>>>> clusters_by_event;
        for (const auto& det_name : target_detectors) {
            auto& by_event = clusters_by_event[det_name];
            for (auto const& cl : batch.clusters.at(det_name)) by_event[cl->eventID()].push_back(cl);
        }
        static const std::vector<std::shared_ptr<Cluster>> empty_clusters;

        size_t n_tracks = batch.tracks.size();
        total_tracks_ += n_tracks;

        size_t n_threads = thread_pool->getThreadCount();
        size_t chunk = (n_tracks + n_threads - 1) / n_threads;
        std::vector<std::future<void>> futures;

        for (size_t i = 0; i < n_tracks; i += chunk) {
            size_t end = std::min(i + chunk, n_tracks);

            futures.push_back(thread_pool->submit([&, i, end]() {
                auto& pm = PlotManager::getInstance();
                size_t l_assoc = 0, l_with_assoc = 0, l_cut_spatial = 0, l_cut_time = 0;
                std::map<std::string, size_t> loc_num_examined;

                for (size_t t_idx = i; t_idx < end; ++t_idx) {
                    auto& track = batch.tracks[t_idx];

                    // Tracked locally (not on Track itself, which only stores
                    // a single closest-cluster pointer with no distance) so
                    // "closest across all target detectors" compares actual
                    // distances consistently.
                    double track_min_dist_sq = std::numeric_limits<double>::max();
                    std::shared_ptr<Cluster> track_closest = nullptr;

                    for (const auto& det_name : target_detectors) {
                        auto& det = geo.getDetector(det_name);
                        std::string dir = getName() + "/" + det_name;
                        auto const& die_labels = die_labels_by_det_.at(det_name);

                        Eigen::Vector3d global_pred = track->interceptAt(det);
                        if (std::isnan(global_pred.x()) || std::isnan(global_pred.y())) continue;
                        Eigen::Vector3d local_pred = det.R_inv * (global_pred - det.position);

                        // Which die the TRACK itself predicts to be on, not
                        // any cluster's - restricts candidates below to
                        // this die only, so a track near the gap between
                        // two dies no longer gets tested against both
                        // sides' boundary-row clusters (each with its own
                        // independent chance to fall inside the ellipse
                        // cut). getColumn()/getRow() take a fractional
                        // position directly - dieOf() rounds to the
                        // nearest pixel internally, no manual rounding
                        // needed here. "" for a detector with no die split
                        // at all (dieOf() returns "" everywhere on it, same
                        // as every cluster's own die - see below), so this
                        // is a no-op there.
                        std::string track_die = det.dieOf(det.getColumn(local_pred.x()), det.getRow(local_pred.y()));

                        int assoc_this_det = 0;
                        std::map<std::string, int> assoc_by_die;
                        std::set<std::string> dies_with_cluster;

                        // Fills a histogram into the combined dir and, if
                        // `die` is non-empty (DetectorGeo::dieOf() of the
                        // cluster this fill belongs to), also into
                        // "<dir>/<die>" - lets a multi-die board's plots be
                        // inspected both combined and per-die without
                        // duplicating every call site below.
                        auto fill1 = [&](const std::string& die, const std::string& name, double v) {
                            pm.fill1D(dir, name, v);
                            if (!die.empty()) pm.fill1D(dir + "/" + die, name, v);
                        };
                        auto fill2 = [&](const std::string& die, const std::string& name, double x, double y) {
                            pm.fill2D(dir, name, x, y);
                            if (!die.empty()) pm.fill2D(dir + "/" + die, name, x, y);
                        };

                        auto& by_event = clusters_by_event.at(det_name);
                        auto ev_it = by_event.find(track->eventID());
                        auto const& det_clusters = (ev_it != by_event.end()) ? ev_it->second : empty_clusters;

                        for (const auto& cluster : det_clusters) {
                            std::string die = det.dieOf(cluster->column(), cluster->row());
                            // Skip a cluster on a different die than the
                            // track's own predicted one - see track_die's
                            // docs above for why. Real dead-zone geometry
                            // (the actual gap between dies) is unaffected:
                            // a track truly in the gap still predicts
                            // SOME die (whichever pixel index it rounds
                            // nearest to), so it's still compared against
                            // that one die's clusters, same as any other
                            // track - it just no longer also gets a second,
                            // independent chance via the other die.
                            if (die != track_die) continue;
                            if (!die.empty()) dies_with_cluster.insert(die);
                            Eigen::Vector3d local_c = det.R_inv * (cluster->global() - det.position);

                            double xdist_centre = std::abs(local_pred.x() - local_c.x());
                            double ydist_centre = std::abs(local_pred.y() - local_c.y());

                            // Distance to the nearest pixel within the cluster
                            // (corry's default association metric).
                            double xdist_nearest = std::numeric_limits<double>::max();
                            double ydist_nearest = std::numeric_limits<double>::max();
                            for (auto const& px : cluster->pixels()) {
                                double lx = det.getLocalX(px->column());
                                double ly = det.getLocalY(px->row());
                                xdist_nearest = std::min(xdist_nearest, std::abs(local_pred.x() - lx));
                                ydist_nearest = std::min(ydist_nearest, std::abs(local_pred.y() - ly));
                            }

                            fill1(die, "hDistXClusterClosestPx", (xdist_centre - xdist_nearest) * 1000.0);
                            fill1(die, "hDistYClusterClosestPx", (ydist_centre - ydist_nearest) * 1000.0);

                            // "Before cuts" 2D distance views, filled for
                            // every examined (track, cluster) pair.
                            fill2(die, "hDist_trackCluster_2D", local_pred.x() - local_c.x(), local_pred.y() - local_c.y());
                            fill2(die, "hDist_trackCluster_2D_huge", local_pred.x() - local_c.x(), local_pred.y() - local_c.y());
                            fill2(die, "hDist_trackCluster_2D_huge_globL", global_pred.x() - cluster->global().x(), global_pred.y() - cluster->global().y());

                            double xdist = use_cluster_centre_ ? xdist_centre : xdist_nearest;
                            double ydist = use_cluster_centre_ ? ydist_centre : ydist_nearest;

                            loc_num_examined[det_name]++;
                            if (!die.empty()) loc_num_examined[det_name + "/" + die]++;

                            // Elliptical cut, matching corry.
                            double norm = (xdist * xdist) / (spatial_cut_x_ * spatial_cut_x_) +
                                          (ydist * ydist) / (spatial_cut_y_ * spatial_cut_y_);
                            if (norm > 1.0) {
                                fill1(die, "hCutHisto", 1);
                                l_cut_spatial++;
                                continue;
                            }

                            if (std::abs(cluster->timestamp() - track->timestamp()) > time_cut_abs_) {
                                fill1(die, "hCutHisto", 2);
                                l_cut_time++;
                                continue;
                            }

                            // Associate every cluster within the cut (not just
                            // the closest one), matching corry.
                            track->addAssociatedCluster(cluster);
                            assoc_this_det++;
                            if (!die.empty()) assoc_by_die[die]++;
                            l_assoc++;

                            fill2(die, "hDist_trackCluster_2D_assoc", local_pred.x() - local_c.x(), local_pred.y() - local_c.y());
                            fill2(die, "hDist_trackCluster_2D_assoc_global", global_pred.x() - cluster->global().x(), global_pred.y() - cluster->global().y());
                            // Matches corry: filled once per associated
                            // cluster (not deduplicated per track) - a track
                            // that associates to several clusters within the
                            // cut fills this once per association, same as
                            // corry's DUTAssociation.cpp does.
                            fill2(die, "hTrackPosAssoc2D", global_pred.x(), global_pred.y());
                            fill2(die, "hClusterPosAssoc2D", local_c.x(), local_c.y());

                            double dist_sq = xdist * xdist + ydist * ydist;
                            if (dist_sq < track_min_dist_sq) {
                                track_min_dist_sq = dist_sq;
                                track_closest = cluster;
                            }
                        }

                        pm.fill1D(dir, "no_assoc_cls", static_cast<double>(assoc_this_det));
                        // Matches corry: hTrackPosNoAssoc2D is only filled
                        // when this event actually had at least one cluster
                        // to attempt association against (a real "found a
                        // cluster nearby but it didn't pass the cuts"), not
                        // for tracks where the event simply had zero
                        // clusters for this detector at all - corry's own
                        // DUTAssociation::run() does `if(clusters.empty())
                        // { hNoAssocCls->Fill(0); continue; }` before ever
                        // reaching this fill.
                        if (assoc_this_det == 0 && !det_clusters.empty()) {
                            pm.fill2D(dir, "hTrackPosNoAssoc2D", global_pred.x(), global_pred.y());
                        }
                        for (auto const& die : die_labels) {
                            auto it = assoc_by_die.find(die);
                            int assoc_n = (it != assoc_by_die.end()) ? it->second : 0;
                            pm.fill1D(dir + "/" + die, "no_assoc_cls", static_cast<double>(assoc_n));
                            if (assoc_n == 0 && dies_with_cluster.count(die)) {
                                pm.fill2D(dir + "/" + die, "hTrackPosNoAssoc2D", global_pred.x(), global_pred.y());
                            }
                        }
                    }

                    if (track_closest) {
                        track->setClosestCluster(track_closest);
                        l_with_assoc++;
                    }
                }

                total_associations_ += l_assoc;
                tracks_with_assoc_ += l_with_assoc;
                cut_spatial_ += l_cut_spatial;
                cut_timing_ += l_cut_time;

                std::lock_guard<std::mutex> lock(counter_mtx_);
                for (auto const& [det_name, count] : loc_num_examined) num_examined_by_det_[det_name] += count;
            }));
        }
        for (auto& f : futures) f.get();
    }

    void DUTAssociation::finalize() {
        // Normalize hCutHisto into a fraction-of-examined-clusters, matching
        // corry's hCutHisto->Scale(1 / double(num_cluster)).
        auto& pm = PlotManager::getInstance();
        for (auto const& [det_name, num_cluster] : num_examined_by_det_) {
            if (num_cluster == 0) continue;
            std::string dir = getName() + "/" + det_name;
            auto& hist = pm.getPlot1D(dir, "hCutHisto");
            std::vector<double> scaled(hist.getBins());
            for (size_t b = 0; b < hist.getBins(); ++b) scaled[b] = hist.getBinContent(b) / static_cast<double>(num_cluster);
            // setData() re-derives entries_ from summing the new bin
            // content (correct for a pure derived/ratio plot with no prior
            // real fills), but here the histogram was already populated by
            // real fill() calls before this rescale - the normalized sum is
            // <= 1.0 so entries_ would truncate to 0. Preserve the real
            // pre-rescale entries count explicitly.
            size_t real_entries = hist.getEntries();
            hist.setData(scaled);
            hist.setEntries(real_entries);
        }

        WR_LOG(STATUS, "In total, " + std::to_string(total_associations_) + " clusters are associated to " +
                            std::to_string(tracks_with_assoc_) + " tracks.");
        WR_LOG(INFO, "Number of tracks with at least one associated cluster: " + std::to_string(tracks_with_assoc_) +
                         " vs total number of tracks: " + std::to_string(total_tracks_));
        if (cut_spatial_ > 0 || cut_timing_ > 0) {
            WR_LOG(INFO, "Discarded by spatial cut: " + std::to_string(cut_spatial_) +
                             ", discarded by timing cut: " + std::to_string(cut_timing_));
        }
    }

    REGISTER_MODULE(DUTAssociation)
}
