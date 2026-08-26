/** @file AnalysisEfficiency.cpp */
#include "AnalysisEfficiency.hpp"
#include "core/ModuleFactory.hpp"
#include "core/detector/GeometryManager.hpp"
#include "core/utils/PlotManager.hpp"
#include "core/utils/Logger.hpp"
#include <algorithm>
#include <cmath>
#include <future>
#include <limits>
#include <set>

namespace framework {

    namespace {
        /// Maps `value` to a 0-based bin index over [lo, hi) split into
        /// `n_bins` equal-width bins, or -1 if out of range.
        int binIndex(double value, double lo, double hi, int n_bins) {
            if (n_bins <= 0 || hi <= lo) return -1;
            int idx = static_cast<int>((value - lo) / (hi - lo) * n_bins);
            if (idx < 0 || idx >= n_bins) return -1;
            return idx;
        }

    }

    AnalysisEfficiency::AnalysisEfficiency(Configuration cfg, Configuration g_cfg, ThreadPool* pool)
        : Module(std::move(cfg), std::move(g_cfg), pool) {
        dut_name_ = config.get<std::string>("name", "");
        type_filter_ = config.get<std::string>("type", "");

        chi2_cut_ = config.get<double>("chi2ndof_cut", 3.0);
        required_detectors_ = config.getArray<std::string>("require_associated_cluster_on");

        masked_pixel_distance_cut_ = config.get<int>("masked_pixel_distance_cut", 1);
        spatial_cut_sensoredge_ = config.get<double>("spatial_cut_sensoredge", 1.0);
        fake_rate_method_ = config.get<std::string>("fake_rate_method", "radius");
        fake_rate_distance_ = config.get<double>("fake_rate_distance", 2.0);
        n_charge_bins_ = config.get<int>("n_charge_bins", 1000);
        charge_histo_range_ = config.get<double>("charge_histo_range", 1000.0);
        sigma_moyal_ = config.get<double>("sigma_moyal", 0.013);
        time_cut_frameedge_ = config.get<double>("time_cut_frameedge", 0.0);

        auto bin_sizes = config.getArray<double>("efficiency_2D_histograms_bin_size");
        bin_size_x_ = bin_sizes.size() >= 1 ? bin_sizes[0] : 0.01;
        bin_size_y_ = bin_sizes.size() >= 2 ? bin_sizes[1] : bin_size_x_;
    }

    void AnalysisEfficiency::initialize() {
        auto& geo = GeometryManager::getInstance();

        // corry binds this module to detectors either via "name" (exactly
        // one) or "type" (auto-instantiated once per matching detector) -
        // replicate the latter the same way DUTAssociation does: one
        // instance internally loops over every matching detector.
        if (!dut_name_.empty()) {
            if (!geo.hasDetector(dut_name_)) {
                throw std::runtime_error("AnalysisEfficiency: detector '" + dut_name_ + "' not found in geometry.");
            }
            target_detectors_.push_back(dut_name_);
        } else if (!type_filter_.empty()) {
            for (const auto& name : geo.getDetectorNames()) {
                if (geo.getDetector(name).type == type_filter_) target_detectors_.push_back(name);
            }
        } else {
            throw std::runtime_error("AnalysisEfficiency: must configure either 'name' or 'type'.");
        }

        auto& pm = PlotManager::getInstance();
        for (const auto& det_name : target_detectors_) {
            auto& det = geo.getDetector(det_name);
            WR_LOG(STATUS, "Initializing Efficiency Analysis for DUT: " + det_name);

            // Builds one full DetAccum + plot set under "<key>" - called
            // once for the combined view (`key == det_name`) and, for a
            // split board, once more per die label under
            // "<det_name>/<die>" (see #die_labels_by_det_/run()). All
            // variants share the exact same geometry-derived binning
            // (`det` itself isn't split, only which associated clusters/
            // pixels feed each accumulator - see runTrackLoop()/
            // runFakeRateLoop()).
            auto setupAccum = [&](const std::string& key) {
                std::string dir = getName() + "/" + key;

                DetAccum acc;
                acc.size_x = det.sizeX();
                acc.size_y = det.sizeY();
                acc.n_bins_x = std::max(1, static_cast<int>(std::round(1.4 * acc.size_x / bin_size_x_)));
                acc.n_bins_y = std::max(1, static_cast<int>(std::round(1.4 * acc.size_y / bin_size_y_)));

                acc.eff_num_global = Eigen::MatrixXd::Zero(acc.n_bins_x, acc.n_bins_y);
                acc.eff_den_global = Eigen::MatrixXd::Zero(acc.n_bins_x, acc.n_bins_y);
                acc.eff_num_local = Eigen::MatrixXd::Zero(acc.n_bins_x, acc.n_bins_y);
                acc.eff_den_local = Eigen::MatrixXd::Zero(acc.n_bins_x, acc.n_bins_y);
                acc.profile_count = Eigen::MatrixXd::Zero(acc.n_bins_x, acc.n_bins_y);
                acc.charge_sum = Eigen::MatrixXd::Zero(acc.n_bins_x, acc.n_bins_y);
                acc.moyal_sum = Eigen::MatrixXd::Zero(acc.n_bins_x, acc.n_bins_y);
                acc.clustersize_sum = Eigen::MatrixXd::Zero(acc.n_bins_x, acc.n_bins_y);
                acc.ncols_sum = Eigen::MatrixXd::Zero(acc.n_bins_x, acc.n_bins_y);
                acc.nrows_sum = Eigen::MatrixXd::Zero(acc.n_bins_x, acc.n_bins_y);
                acc.toa_sum = Eigen::MatrixXd::Zero(acc.n_bins_x, acc.n_bins_y);
                acc.fake_pixel_map = Eigen::MatrixXd::Zero(det.n_pixels_x, det.n_pixels_y);
                acc.fake_vs_time_sum.assign(3000, 0.0);
                acc.fake_vs_time_count.assign(3000, 0.0);
                acc.fake_vs_time_long_sum.assign(3000, 0.0);
                acc.fake_vs_time_long_count.assign(3000, 0.0);

                double abs_sx = std::abs(acc.size_x), abs_sy = std::abs(acc.size_y);
                double gx_lo = det.position.x() - 0.7 * abs_sx, gx_hi = det.position.x() + 0.7 * abs_sx;
                double gy_lo = det.position.y() - 0.7 * abs_sy, gy_hi = det.position.y() + 0.7 * abs_sy;
                pm.registerPlot2D(dir, "globalEfficiencyMap_trackPos", acc.n_bins_x, gx_lo, gx_hi, acc.n_bins_y, gy_lo, gy_hi);
                pm.registerPlot2D(dir, "localEfficiencyMap_trackPos", acc.n_bins_x, -0.7 * abs_sx, 0.7 * abs_sx, acc.n_bins_y, -0.7 * abs_sy, 0.7 * abs_sy);
                pm.registerPlot2D(dir, "localChargeMapMean_trackPos", acc.n_bins_x, -0.7 * abs_sx, 0.7 * abs_sx, acc.n_bins_y, -0.7 * abs_sy, 0.7 * abs_sy);
                pm.registerPlot2D(dir, "localMoyalCharge_trackPos", acc.n_bins_x, -0.7 * abs_sx, 0.7 * abs_sx, acc.n_bins_y, -0.7 * abs_sy, 0.7 * abs_sy);
                pm.registerPlot2D(dir, "localMPVCharge_trackPos", acc.n_bins_x, -0.7 * abs_sx, 0.7 * abs_sx, acc.n_bins_y, -0.7 * abs_sy, 0.7 * abs_sy);
                pm.registerPlot2D(dir, "clusterSize_trackPos", acc.n_bins_x, -0.7 * abs_sx, 0.7 * abs_sx, acc.n_bins_y, -0.7 * abs_sy, 0.7 * abs_sy);
                pm.registerPlot2D(dir, "clusterNCols_trackPos", acc.n_bins_x, -0.7 * abs_sx, 0.7 * abs_sx, acc.n_bins_y, -0.7 * abs_sy, 0.7 * abs_sy);
                pm.registerPlot2D(dir, "clusterNRows_trackPos", acc.n_bins_x, -0.7 * abs_sx, 0.7 * abs_sx, acc.n_bins_y, -0.7 * abs_sy, 0.7 * abs_sy);
                pm.registerPlot2D(dir, "ToAMap_trackPos", acc.n_bins_x, -0.7 * abs_sx, 0.7 * abs_sx, acc.n_bins_y, -0.7 * abs_sy, 0.7 * abs_sy);

                double dist_range = std::sqrt(det.pitch_x * det.pitch_y);
                int dist_bins = std::max(1, static_cast<int>(dist_range));
                pm.registerPlot1D(dir, "distanceTrackHit", dist_bins, 0, dist_range);
                pm.registerPlot2D(dir, "distanceTrackHit2D", 150, -1.5 * det.pitch_x, 1.5 * det.pitch_x, 150, -1.5 * det.pitch_y, 1.5 * det.pitch_y);

                pm.registerGraph(dir, "efficiency_vs_event");

                // Fake rate plots (corry's createFakeRatePlots(), unconditional).
                pm.registerPlot1D(dir, "hFakePixelPerEvent", 25, -0.5, 24.5);
                pm.registerPlot2D(dir, "fakePixelPerEventMap", det.n_pixels_x, -0.5, det.n_pixels_x - 0.5, det.n_pixels_y, -0.5, det.n_pixels_y - 0.5);
                pm.registerPlot1D(dir, "fakePixelPerEventVsTime", 3000, 0, 3000);
                // Matches corry's own name for this histogram,
                // "efficiencyVsTimeLong", for output parity.
                pm.registerPlot1D(dir, "efficiencyVsTimeLong", 3000, 0, 30000);
                pm.registerPlot1D(dir, "hFakePixelCharge", n_charge_bins_, 0.0, charge_histo_range_);
                pm.registerPlot1D(dir, "hFakeClusterCharge", n_charge_bins_, 0.0, charge_histo_range_);
                pm.registerPlot1D(dir, "hFakeClusterPerEvent", 25, -0.5, 24.5);
                pm.registerPlot1D(dir, "hFakeClusterSize", 25, -0.5, 24.5);

                accum_[key] = std::move(acc);
            };

            setupAccum(det_name);

            // Uses DetectorGeo::dieLabels() rather than a bare n_pixels_y
            // check: some single-die 2D-grid sensors also satisfy
            // n_pixels_y >= 2, and would otherwise get empty, never-filled
            // per-die plots registered (dieOf(), used below in run(),
            // always returns "" for them - see dieLabels()'s own docs).
            auto labels = det.dieLabels();
            if (!labels.empty()) {
                for (auto const& die : labels) setupAccum(det_name + "/" + die);
                die_labels_by_det_[det_name] = std::move(labels);
            } else if (det.hasUnrecognizedLayout()) {
                WR_LOG(WARNING, "AnalysisEfficiency: detector '" + det_name + "' has unrecognized layout '" +
                                     det.layout + "' - not splitting into per-die plots. Add a dedicated Layout "
                                     "subclass (src/core/detector/layouts/) if this layout needs die-level splitting.");
            }
        }
    }

    void AnalysisEfficiency::run(DataBatch& batch) {
        if (batch.tracks.empty()) return;

        // Built once and shared across every target detector's fake-rate
        // pass below - grouping is independent of detector, so rebuilding
        // it per detector when several share one config block (type = "...")
        // would be redundant.
        std::map<uint64_t, std::vector<std::shared_ptr<Track>>> tracks_by_event;
        for (auto const& t : batch.tracks) tracks_by_event[t->eventID()].push_back(t);

        for (const auto& det_name : target_detectors_) {
            runTrackLoop(batch, det_name, det_name, "");
            runFakeRateLoop(batch, det_name, det_name, "", tracks_by_event);

            auto it = die_labels_by_det_.find(det_name);
            if (it != die_labels_by_det_.end()) {
                for (auto const& die : it->second) {
                    std::string accum_key = det_name + "/" + die;
                    runTrackLoop(batch, det_name, accum_key, die);
                    runFakeRateLoop(batch, det_name, accum_key, die, tracks_by_event);
                }
            }
        }
    }

    void AnalysisEfficiency::runTrackLoop(DataBatch& batch, const std::string& det_name, const std::string& accum_key,
                                           const std::string& die_filter) {
        auto& geo = GeometryManager::getInstance();
        auto& det = geo.getDetector(det_name);
        std::string dir = getName() + "/" + accum_key;

        int n_bins_x, n_bins_y;
        double abs_sx, abs_sy;
        {
            auto& acc = accum_.at(accum_key);
            n_bins_x = acc.n_bins_x; n_bins_y = acc.n_bins_y;
            abs_sx = std::abs(acc.size_x); abs_sy = std::abs(acc.size_y);
        }

        size_t n_tracks = batch.tracks.size();
        size_t n_threads = thread_pool->getThreadCount();
        size_t chunk = (n_tracks + n_threads - 1) / n_threads;
        std::vector<std::future<void>> futures;

        for (size_t i = 0; i < n_tracks; i += chunk) {
            size_t end = std::min(i + chunk, n_tracks);
            futures.push_back(thread_pool->submit([&, i, end]() {
                auto& pm = PlotManager::getInstance();

                Eigen::MatrixXd loc_num_g = Eigen::MatrixXd::Zero(n_bins_x, n_bins_y), loc_den_g = Eigen::MatrixXd::Zero(n_bins_x, n_bins_y);
                Eigen::MatrixXd loc_num_l = Eigen::MatrixXd::Zero(n_bins_x, n_bins_y), loc_den_l = Eigen::MatrixXd::Zero(n_bins_x, n_bins_y);
                Eigen::MatrixXd loc_prof_count = Eigen::MatrixXd::Zero(n_bins_x, n_bins_y);
                Eigen::MatrixXd loc_charge = Eigen::MatrixXd::Zero(n_bins_x, n_bins_y);
                Eigen::MatrixXd loc_moyal = Eigen::MatrixXd::Zero(n_bins_x, n_bins_y);
                Eigen::MatrixXd loc_clsize = Eigen::MatrixXd::Zero(n_bins_x, n_bins_y);
                Eigen::MatrixXd loc_ncols = Eigen::MatrixXd::Zero(n_bins_x, n_bins_y);
                Eigen::MatrixXd loc_nrows = Eigen::MatrixXd::Zero(n_bins_x, n_bins_y);
                Eigen::MatrixXd loc_toa = Eigen::MatrixXd::Zero(n_bins_x, n_bins_y);

                size_t l_track = 0, l_dut = 0, l_chi2 = 0, l_masked = 0, l_requirecluster = 0, l_total = 0, l_matched = 0;
                std::map<uint64_t, std::pair<size_t, size_t>> loc_event_counts;

                for (size_t t_idx = i; t_idx < end; ++t_idx) {
                    auto const& track = batch.tracks[t_idx];
                    l_track++;

                    Eigen::Vector3d global_pred = track->interceptAt(det);
                    if (std::isnan(global_pred.x()) || std::isnan(global_pred.y())) { l_dut++; continue; }
                    Eigen::Vector3d local_pred = det.R_inv * (global_pred - det.position);

                    if (!det.hasIntercept(local_pred.x(), local_pred.y(), spatial_cut_sensoredge_)) { l_dut++; continue; }
                    // For a per-die accumulator, a track whose own intercept
                    // lands on the other die still passes the whole-detector
                    // hasIntercept() check above (the combined footprint,
                    // both dies together), but can never match this die's
                    // cluster filter below. Treat it as outside the DUT here
                    // (same bucket as the check above) so it doesn't inflate
                    // this die's denominator - left unhandled, per-die
                    // efficiency would be capped near 50%. The combined
                    // (die_filter-empty) pass is unaffected since it never
                    // applies this check.
                    if (!die_filter.empty() &&
                        det.dieOf(det.getColumn(local_pred.x()), det.getRow(local_pred.y())) != die_filter) {
                        l_dut++;
                        continue;
                    }
                    if (track->getChi2ndof() > chi2_cut_) { l_chi2++; continue; }
                    // isIsolated(): always true for now, see Track::isIsolated docs.
                    if (!track->isIsolated()) continue;
                    // isWithinROI: no `continue` in corry either (would only
                    // gate the in-pixel-ROI plots, which aren't reproduced
                    // here) - nothing further to do.
                    if (det.hitMasked(local_pred.x(), local_pred.y(), masked_pixel_distance_cut_)) { l_masked++; continue; }
                    // frame-edge cut: no-op, see time_cut_frameedge_ docs.

                    bool requirements_met = true;
                    for (const auto& req : required_detectors_) {
                        bool has_req = false;
                        for (auto const& c : track->getAssociatedClusters()) {
                            if (c->detectorID() == req) { has_req = true; break; }
                        }
                        if (!has_req) { requirements_met = false; break; }
                    }
                    if (!requirements_met) { l_requirecluster++; continue; }

                    l_total++;

                    bool has_associated = false;
                    for (auto const& c : track->getAssociatedClusters()) {
                        if (c->detectorID() != det_name) continue;
                        if (!die_filter.empty() && det.dieOf(c->column(), c->row()) != die_filter) continue;
                        has_associated = true;
                        break;
                    }
                    if (has_associated) l_matched++;

                    int gcol = binIndex(global_pred.x(), det.position.x() - 0.7 * abs_sx, det.position.x() + 0.7 * abs_sx, n_bins_x);
                    int grow = binIndex(global_pred.y(), det.position.y() - 0.7 * abs_sy, det.position.y() + 0.7 * abs_sy, n_bins_y);
                    int lcol = binIndex(local_pred.x(), -0.7 * abs_sx, 0.7 * abs_sx, n_bins_x);
                    int lrow = binIndex(local_pred.y(), -0.7 * abs_sy, 0.7 * abs_sy, n_bins_y);

                    if (gcol >= 0 && grow >= 0) { loc_den_g(gcol, grow) += 1; if (has_associated) loc_num_g(gcol, grow) += 1; }
                    if (lcol >= 0 && lrow >= 0) { loc_den_l(lcol, lrow) += 1; if (has_associated) loc_num_l(lcol, lrow) += 1; }

                    auto& ev_count = loc_event_counts[track->eventID()];
                    ev_count.first++;
                    if (has_associated) ev_count.second++;

                    if (has_associated) {
                        // Closest associated cluster ON THIS DETECTOR (found
                        // locally rather than via Track::getClosestCluster(),
                        // which has no per-detector filter and could belong
                        // to a different detector when several DUTs of the
                        // same type share one DUTAssociation instance).
                        std::shared_ptr<Cluster> closest;
                        double closest_d2 = std::numeric_limits<double>::max();
                        for (auto const& c : track->getAssociatedClusters()) {
                            if (c->detectorID() != det_name) continue;
                            if (!die_filter.empty() && det.dieOf(c->column(), c->row()) != die_filter) continue;
                            Eigen::Vector3d local_c = det.R_inv * (c->global() - det.position);
                            double dx = local_pred.x() - local_c.x(), dy = local_pred.y() - local_c.y();
                            double d2 = dx * dx + dy * dy;
                            if (d2 < closest_d2) { closest_d2 = d2; closest = c; }
                        }
                        if (closest) {
                            Eigen::Vector3d local_c = det.R_inv * (closest->global() - det.position);
                            double dx = local_pred.x() - local_c.x(), dy = local_pred.y() - local_c.y();
                            pm.fill2D(dir, "distanceTrackHit2D", dx, dy);
                            pm.fill1D(dir, "distanceTrackHit", std::sqrt(dx * dx + dy * dy));
                        }

                        for (auto const& c : track->getAssociatedClusters()) {
                            if (c->detectorID() != det_name) continue;
                            if (!die_filter.empty() && det.dieOf(c->column(), c->row()) != die_filter) continue;
                            if (lcol < 0 || lrow < 0) continue;
                            loc_prof_count(lcol, lrow) += 1;
                            loc_charge(lcol, lrow) += c->charge();
                            loc_moyal(lcol, lrow) += std::exp(-c->charge() / sigma_moyal_);
                            loc_clsize(lcol, lrow) += static_cast<double>(c->size());
                            loc_ncols(lcol, lrow) += static_cast<double>(c->columnWidth());
                            loc_nrows(lcol, lrow) += static_cast<double>(c->rowWidth());
                            loc_toa(lcol, lrow) += c->timestamp();
                        }
                    }
                }

                std::lock_guard<std::mutex> lock(acc_mtx_);
                auto& shared = accum_.at(accum_key);
                shared.eff_num_global += loc_num_g; shared.eff_den_global += loc_den_g;
                shared.eff_num_local += loc_num_l; shared.eff_den_local += loc_den_l;
                shared.profile_count += loc_prof_count;
                shared.charge_sum += loc_charge; shared.moyal_sum += loc_moyal;
                shared.clustersize_sum += loc_clsize; shared.ncols_sum += loc_ncols; shared.nrows_sum += loc_nrows;
                shared.toa_sum += loc_toa;
                shared.n_track += l_track; shared.n_dut += l_dut; shared.n_chi2 += l_chi2;
                shared.n_masked += l_masked; shared.n_requirecluster += l_requirecluster;
                shared.total_tracks += l_total; shared.matched_tracks += l_matched;
                for (auto const& [ev, cnt] : loc_event_counts) {
                    shared.per_event_counts[ev].first += cnt.first;
                    shared.per_event_counts[ev].second += cnt.second;
                }
            }));
        }
        for (auto& f : futures) f.get();
    }

    void AnalysisEfficiency::runFakeRateLoop(DataBatch& batch, const std::string& det_name, const std::string& accum_key,
                                              const std::string& die_filter,
                                              const std::map<uint64_t, std::vector<std::shared_ptr<Track>>>& tracks_by_event) {
        auto& geo = GeometryManager::getInstance();
        auto& det = geo.getDetector(det_name);
        std::string dir = getName() + "/" + accum_key;

        std::map<uint64_t, std::vector<std::shared_ptr<Cluster>>> clusters_by_event;
        if (batch.clusters.count(det_name)) {
            for (auto const& c : batch.clusters.at(det_name)) {
                if (!die_filter.empty() && det.dieOf(c->column(), c->row()) != die_filter) continue;
                clusters_by_event[c->eventID()].push_back(c);
            }
        }
        std::map<uint64_t, std::vector<std::shared_ptr<Pixel>>> pixels_by_event;
        if (fake_rate_method_ == "edge" && batch.pixels.count(det_name)) {
            for (auto const& p : batch.pixels.at(det_name)) {
                if (!die_filter.empty() && det.dieOf(p->column(), p->row()) != die_filter) continue;
                pixels_by_event[p->eventID()].push_back(p);
            }
        }

        std::set<uint64_t> event_id_set;
        for (auto const& [ev, v] : tracks_by_event) event_id_set.insert(ev);
        for (auto const& [ev, v] : clusters_by_event) event_id_set.insert(ev);
        std::vector<uint64_t> event_list(event_id_set.begin(), event_id_set.end());

        size_t n_events = event_list.size();
        size_t n_threads = thread_pool->getThreadCount();
        size_t chunk = (n_events + n_threads - 1) / n_threads;
        std::vector<std::future<void>> futures;
        static const std::vector<std::shared_ptr<Track>> empty_tracks;
        static const std::vector<std::shared_ptr<Cluster>> empty_clusters;
        static const std::vector<std::shared_ptr<Pixel>> empty_pixels;

        for (size_t i = 0; i < n_events; i += chunk) {
            size_t end = std::min(i + chunk, n_events);
            futures.push_back(thread_pool->submit([&, i, end]() {
                auto& pm = PlotManager::getInstance();
                Eigen::MatrixXd loc_fake_map = Eigen::MatrixXd::Zero(det.n_pixels_x, det.n_pixels_y);
                std::vector<double> loc_t_sum(3000, 0.0), loc_t_count(3000, 0.0);
                std::vector<double> loc_tl_sum(3000, 0.0), loc_tl_count(3000, 0.0);
                size_t loc_events = 0;

                for (size_t idx = i; idx < end; ++idx) {
                    uint64_t ev = event_list[idx];
                    auto const& tracks = tracks_by_event.count(ev) ? tracks_by_event.at(ev) : empty_tracks;
                    auto const& clusters = clusters_by_event.count(ev) ? clusters_by_event.at(ev) : empty_clusters;

                    // Event "start time" proxy (s): Warlock has no
                    // Event/frame-boundary object, so the earliest track
                    // timestamp in this event stands in for corry's
                    // event->start(). Native timestamps are ns (Units.cpp).
                    double event_time_s = 0.0;
                    bool have_time = false;
                    for (auto const& t : tracks) {
                        if (!have_time || t->timestamp() < event_time_s) { event_time_s = t->timestamp(); have_time = true; }
                    }
                    event_time_s /= 1e9;

                    int fake_hits = 0, fake_clusters = 0;
                    bool do_fill_summary = true;

                    if (fake_rate_method_ == "edge") {
                        bool track_in_active = false;
                        for (auto const& t : tracks) {
                            Eigen::Vector3d global_pred = t->interceptAt(det);
                            if (std::isnan(global_pred.x())) continue;
                            Eigen::Vector3d local_pred = det.R_inv * (global_pred - det.position);
                            if (det.hasIntercept(local_pred.x(), local_pred.y(), -fake_rate_distance_)) { track_in_active = true; break; }
                        }
                        if (track_in_active) {
                            do_fill_summary = false;
                        } else {
                            auto const& pixels = pixels_by_event.count(ev) ? pixels_by_event.at(ev) : empty_pixels;
                            for (auto const& px : pixels) {
                                fake_hits++;
                                pm.fill1D(dir, "hFakePixelCharge", px->charge());
                                if (px->column() >= 0 && px->column() < det.n_pixels_x && px->row() >= 0 && px->row() < det.n_pixels_y) {
                                    loc_fake_map(px->column(), px->row()) += 1;
                                }
                            }
                            for (auto const& cl : clusters) {
                                fake_clusters++;
                                pm.fill1D(dir, "hFakeClusterCharge", cl->charge());
                                pm.fill1D(dir, "hFakeClusterSize", static_cast<double>(cl->size()));
                            }
                        }
                    } else { // radius (default)
                        for (auto const& cl : clusters) {
                            Eigen::Vector3d local_c = det.R_inv * (cl->global() - det.position);
                            bool track_too_close = false;
                            for (auto const& t : tracks) {
                                Eigen::Vector3d global_pred = t->interceptAt(det);
                                if (std::isnan(global_pred.x())) continue;
                                Eigen::Vector3d local_pred = det.R_inv * (global_pred - det.position);
                                // Ported literally from corry: tracks whose
                                // (widened-by-fake_rate_distance) intercept
                                // still lies within bounds are skipped from
                                // this proximity check entirely.
                                if (det.hasIntercept(local_pred.x(), local_pred.y(), -fake_rate_distance_)) continue;

                                double norm_x = (local_pred.x() - local_c.x()) / det.pitch_x;
                                double norm_y = (local_pred.y() - local_c.y()) / det.pitch_y;
                                if (std::sqrt(norm_x * norm_x + norm_y * norm_y) < fake_rate_distance_) { track_too_close = true; break; }
                            }
                            if (!track_too_close) {
                                fake_clusters++;
                                pm.fill1D(dir, "hFakeClusterCharge", cl->charge());
                                pm.fill1D(dir, "hFakeClusterSize", static_cast<double>(cl->size()));
                                for (auto const& px : cl->pixels()) {
                                    fake_hits++;
                                    pm.fill1D(dir, "hFakePixelCharge", px->charge());
                                    if (px->column() >= 0 && px->column() < det.n_pixels_x && px->row() >= 0 && px->row() < det.n_pixels_y) {
                                        loc_fake_map(px->column(), px->row()) += 1;
                                    }
                                }
                            }
                        }
                    }

                    if (do_fill_summary) {
                        pm.fill1D(dir, "hFakePixelPerEvent", static_cast<double>(fake_hits));
                        pm.fill1D(dir, "hFakeClusterPerEvent", static_cast<double>(fake_clusters));
                        loc_events++;

                        if (have_time) {
                            int tbin = static_cast<int>(event_time_s);
                            if (tbin >= 0 && tbin < 3000) { loc_t_sum[tbin] += fake_hits; loc_t_count[tbin] += 1; }
                            int tbin_long = static_cast<int>(event_time_s / 10.0);
                            if (tbin_long >= 0 && tbin_long < 3000) { loc_tl_sum[tbin_long] += fake_hits; loc_tl_count[tbin_long] += 1; }
                        }
                    }
                }

                std::lock_guard<std::mutex> lock(acc_mtx_);
                auto& shared = accum_.at(accum_key);
                shared.fake_pixel_map += loc_fake_map;
                shared.total_events += loc_events;
                for (int b = 0; b < 3000; ++b) {
                    shared.fake_vs_time_sum[b] += loc_t_sum[b]; shared.fake_vs_time_count[b] += loc_t_count[b];
                    shared.fake_vs_time_long_sum[b] += loc_tl_sum[b]; shared.fake_vs_time_long_count[b] += loc_tl_count[b];
                }
            }));
        }
        for (auto& f : futures) f.get();
    }

    void AnalysisEfficiency::finalize() {
        auto& pm = PlotManager::getInstance();

        // Turns one DetAccum's raw sums into its final ratio/mean plots and
        // logs its track-selection flow - shared by the combined
        // (`key == det_name`) accumulator and, for split boards, its
        // per-die siblings.
        auto fillAccumPlots = [&](const std::string& key) {
            auto& acc = accum_.at(key);
            std::string dir = getName() + "/" + key;

            auto ratio = [](const Eigen::MatrixXd& num, const Eigen::MatrixXd& den) {
                return num.binaryExpr(den, [](double n, double d) { return d > 0 ? n / d : 0.0; });
            };
            pm.getPlot2D(dir, "globalEfficiencyMap_trackPos").setData(ratio(acc.eff_num_global, acc.eff_den_global));
            pm.getPlot2D(dir, "localEfficiencyMap_trackPos").setData(ratio(acc.eff_num_local, acc.eff_den_local));
            pm.getPlot2D(dir, "clusterSize_trackPos").setData(ratio(acc.clustersize_sum, acc.profile_count));
            pm.getPlot2D(dir, "clusterNCols_trackPos").setData(ratio(acc.ncols_sum, acc.profile_count));
            pm.getPlot2D(dir, "clusterNRows_trackPos").setData(ratio(acc.nrows_sum, acc.profile_count));
            pm.getPlot2D(dir, "ToAMap_trackPos").setData(ratio(acc.toa_sum, acc.profile_count));
            pm.getPlot2D(dir, "localChargeMapMean_trackPos").setData(ratio(acc.charge_sum, acc.profile_count));

            // Recompute the Moyal charge map into its MPV (matches corry's
            // finalize(): both localMoyalCharge_trackPos AND
            // localMPVCharge_trackPos end up holding the same recomputed
            // charge value, -sigma*log(mean(exp(-Q/sigma)))).
            Eigen::MatrixXd moyal_mean = ratio(acc.moyal_sum, acc.profile_count);
            Eigen::MatrixXd mpv_charge = Eigen::MatrixXd::Zero(acc.n_bins_x, acc.n_bins_y);
            for (int ix = 0; ix < acc.n_bins_x; ++ix) {
                for (int iy = 0; iy < acc.n_bins_y; ++iy) {
                    double m = moyal_mean(ix, iy);
                    if (m > 0) mpv_charge(ix, iy) = -sigma_moyal_ * std::log(m);
                }
            }
            pm.getPlot2D(dir, "localMoyalCharge_trackPos").setData(mpv_charge);
            pm.getPlot2D(dir, "localMPVCharge_trackPos").setData(mpv_charge);

            // Fake pixel map, normalized by number of events actually
            // contributing an hFakePixelPerEvent entry (corry:
            // fakePixelPerEventMap->Scale(1./hFakePixelPerEvent->GetEntries())).
            if (acc.total_events > 0) {
                Eigen::MatrixXd fake_map_scaled = acc.fake_pixel_map / static_cast<double>(acc.total_events);
                pm.getPlot2D(dir, "fakePixelPerEventMap").setData(fake_map_scaled);
            }

            auto& t_hist = pm.getPlot1D(dir, "fakePixelPerEventVsTime");
            std::vector<double> t_mean(3000, 0.0);
            for (int b = 0; b < 3000; ++b) if (acc.fake_vs_time_count[b] > 0) t_mean[b] = acc.fake_vs_time_sum[b] / acc.fake_vs_time_count[b];
            t_hist.setData(t_mean);

            auto& tl_hist = pm.getPlot1D(dir, "efficiencyVsTimeLong");
            std::vector<double> tl_mean(3000, 0.0);
            for (int b = 0; b < 3000; ++b) if (acc.fake_vs_time_long_count[b] > 0) tl_mean[b] = acc.fake_vs_time_long_sum[b] / acc.fake_vs_time_long_count[b];
            tl_hist.setData(tl_mean);

            // Efficiency-vs-event graph: one point per event (rather than
            // corry's one point per accepted track within the event, which
            // isn't reconstructible from batched/multi-threaded processing)
            // holding the cumulative efficiency up to and including that
            // event, in event order.
            size_t cum_total = 0, cum_matched = 0;
            for (auto const& [ev, counts] : acc.per_event_counts) {
                cum_total += counts.first;
                cum_matched += counts.second;
                double eff = cum_total > 0 ? static_cast<double>(cum_matched) / static_cast<double>(cum_total) : 0.0;
                pm.addPoint(dir, "efficiency_vs_event", static_cast<double>(ev), eff);
            }

            double total_eff = acc.total_tracks > 0 ? 100.0 * static_cast<double>(acc.matched_tracks) / static_cast<double>(acc.total_tracks) : 0.0;

            WR_LOG(STATUS, "Track selection flow for " + key + ": " + std::to_string(acc.n_track) +
                                " total, -" + std::to_string(acc.n_dut) + " outside DUT, -" + std::to_string(acc.n_chi2) +
                                " by chi2, -" + std::to_string(acc.n_masked) + " near masked pixel, -" +
                                std::to_string(acc.n_requirecluster) + " missing required cluster, " +
                                std::to_string(acc.total_tracks) + " accepted.");
            WR_LOG(STATUS, "Total efficiency of detector " + key + ": " + std::to_string(total_eff) + "%, measured with " +
                               std::to_string(acc.matched_tracks) + "/" + std::to_string(acc.total_tracks) + " matched/total tracks");
        };

        for (const auto& det_name : target_detectors_) {
            fillAccumPlots(det_name);
            auto it = die_labels_by_det_.find(det_name);
            if (it != die_labels_by_det_.end()) {
                for (auto const& die : it->second) fillAccumPlots(det_name + "/" + die);
            }
        }
    }

    REGISTER_MODULE(AnalysisEfficiency)
}
