/** @file AnalysisTiming.cpp */
#include "AnalysisTiming.hpp"
#include "core/ModuleFactory.hpp"
#include "core/detector/GeometryManager.hpp"
#include "core/utils/PlotManager.hpp"
#include "core/utils/Logger.hpp"
#include "objects/Cluster.hpp"
#include <future>
#include <cmath>

namespace framework {
    /// Shared by run()'s live fill path and finalize()'s deferred (time-walk
    /// fit) replay path, so the two can never label the same die pair
    /// differently.
    static std::string dieLabel(const std::string& name, const std::string& die) {
        return die.empty() ? name : name + "_" + die;
    }

    AnalysisTiming::AnalysisTiming(Configuration cfg, Configuration g_cfg, ThreadPool* pool)
        : Module(std::move(cfg), std::move(g_cfg), pool) {
        
        det_names_ = config.getArray<std::string>("detectors");
        required_detectors_ = config.getArray<std::string>("require_associated_cluster_on");

        // Methods: "none", "inverse_charge", "inverse_sqrt_charge"
        time_walk_method_ = config.get<std::string>("time_walk_method", "none");
        require_same_trigger_group_ = config.get<bool>("require_same_trigger_group", true);
        for (const auto& det : det_names_) {
            time_walk_a_[det] = config.get<double>("time_walk_a_" + det, 0.0);
            time_walk_b_[det] = config.get<double>("time_walk_b_" + det, 0.0);
        }

        time_walk_correction_ = config.get<bool>("time_walk_correction", false);
        if (time_walk_correction_ && time_walk_method_ == "none") {
            WR_LOG(WARNING, "AnalysisTiming: time_walk_correction=true but time_walk_method is \"none\" - "
                             "there's no functional form to fit against, so this has nothing to do. Set "
                             "time_walk_method to \"inverse_charge\" or \"inverse_sqrt_charge\" too. Disabling "
                             "time_walk_correction for this run.");
            time_walk_correction_ = false;
        }
    }

    void AnalysisTiming::initialize() {
        auto& pm = PlotManager::getInstance();
        auto& geo = GeometryManager::getInstance();
        WR_LOG(STATUS, "Initializing AnalysisTiming for " + std::to_string(det_names_.size()) + " detectors.");

        // Plot names match corry's own AnalysisTiming module exactly, so
        // they're directly comparable entry-for-entry. Binning ranges do
        // NOT copy corry's literal initial values (charge [0,100], time
        // [-1,1]): those are ROOT TH1's starting point before
        // SetCanExtend() auto-widens them as data arrives, a feature
        // Warlock's fixed-range histograms don't have. Charge and delta-ToA
        // scales here instead use a wide fixed range up front, since a
        // detector pair on a different clock reference can span a delta-ToA
        // scale an order of magnitude wider than pairs sharing one - a
        // range too narrow for that would silently drop those fills while
        // still counting them as entries (Histogram1D::fill's entries_++
        // happens unconditionally, matching ROOT's own out-of-range Fill()
        // behavior, but bin content only updates in-range). +/-5 comfortably
        // covers every pair with ~1ps bins either way.

        // Die split, applied uniformly to any real multi-die board
        // (DetectorGeo::dieLabels()), same convention as
        // DUTAssociation/AnalysisEfficiency - a detector without a
        // recognized multi-die layout (or not even a pixel detector, e.g.
        // a reference plane) only ever gets the "" (unsplit) variant, so
        // its plot names are completely unchanged from before the split.
        // run() fills every "" (combined) plot exactly as before AND, for
        // a cluster whose own die is known, the matching "_<die>" variant
        // registered here - purely additive.
        auto dieVariants = [&](const std::string& name) -> std::vector<std::string> {
            if (!geo.hasDetector(name)) return {""};
            const auto& det = geo.getDetector(name);
            // Uses DetectorGeo::dieLabels() rather than a bare n_pixels_y
            // check: some single-die 2D-grid sensors also satisfy
            // n_pixels_y >= 2, and would otherwise get empty, never-filled
            // per-die plots registered (dieOf() always returns "" for them
            // - see dieLabels()'s own docs).
            auto labels = det.dieLabels();
            if (!labels.empty()) {
                labels.insert(labels.begin(), "");
                return labels;
            }
            if (det.hasUnrecognizedLayout()) {
                WR_LOG(WARNING, "AnalysisTiming: detector '" + name + "' has unrecognized layout '" +
                                     det.layout + "' - not splitting into per-die plots. Add a dedicated Layout "
                                     "subclass (src/core/detector/layouts/) if this layout needs die-level splitting.");
            }
            return {""};
        };
        for (size_t i = 0; i < det_names_.size(); ++i) {
            const std::string& det1 = det_names_[i];
            auto dies1 = dieVariants(det1);

            for (auto const& d1 : dies1) {
                std::string label1 = dieLabel(det1, d1);
                // Y-axis is raw earliestPixelTimestamp() (event baseline
                // time + a small CFD offset, see
                // Cluster::earliestPixelTimestamp()) - an absolute value
                // that grows across the whole run, unlike the small,
                // relative few-ns scale dtoa_*/timewalk_* use below (those
                // are differences between two simultaneous clusters, this
                // is one cluster's raw, un-differenced timestamp). 1e13
                // gives a comfortable margin over a long run; resolution is
                // necessarily coarse at this scale, since this plot is a
                // rough "does charge correlate with raw arrival time"
                // sanity check, not a precision tool - timewalk_*_in_* is
                // the right plot for that.
                pm.registerPlot2D(getName(), "charge_toa_" + label1, 500, 0, 1.2, 1000, 0.0, 1e13);
                pm.registerPlot1D(getName(), "npairs_" + label1, 10, -0.5, 9.5);
                pm.registerPlot1D(getName(), "ncat_" + label1, 10, -0.5, 9.5);
            }

            for (size_t j = i + 1; j < det_names_.size(); ++j) {
                const std::string& det2 = det_names_[j];

                // Always-filled combined pair (both sides unqualified,
                // pools every die on both sides) - registered once,
                // unconditionally, matching run()'s own unconditional fill.
                std::string combined_pair = det1 + "_vs_" + det2;
                pm.registerPlot1D(getName(), "dtoa_" + combined_pair, 10000, -5.0, 5.0);
                pm.registerPlot2D(getName(), "timewalk_" + det1 + "_in_" + combined_pair, 100, 0, 1.2, 100, -5.0, 5.0);

                // Fully-qualified per-die variants only - deliberately NOT
                // the full dies1 x dies2 cross product. A cluster's own die
                // (DetectorGeo::dieOf()) is never "" once that detector
                // splits, so run() can never fill a "one side qualified,
                // other left blank" combination for a split detector -
                // registering it anyway would just create a
                // permanently-empty, misleadingly-named plot. If a detector
                // doesn't split at all, its own side simply stays "" (not
                // ambiguous - there's only one die to mean).
                std::vector<std::string> reachable1 = geo.hasDetector(det1) ? geo.getDetector(det1).dieLabels() : std::vector<std::string>{};
                std::vector<std::string> reachable2 = geo.hasDetector(det2) ? geo.getDetector(det2).dieLabels() : std::vector<std::string>{};
                if (reachable1.empty()) reachable1 = {""};
                if (reachable2.empty()) reachable2 = {""};
                for (auto const& d1 : reachable1) {
                    for (auto const& d2 : reachable2) {
                        if (d1.empty() && d2.empty()) continue; // already registered as combined_pair above
                        std::string label1 = dieLabel(det1, d1);
                        std::string label2 = dieLabel(det2, d2);
                        std::string pair_variant = label1 + "_vs_" + label2;

                        pm.registerPlot1D(getName(), "dtoa_" + pair_variant, 10000, -5.0, 5.0);

                        // Warlock-only addition, not in corry: amplitude vs
                        // delta-ToA correlation, useful for deriving time-walk
                        // parameters offline without corry's per-pair TTree.
                        pm.registerPlot2D(getName(), "timewalk_" + label1 + "_in_" + pair_variant, 100, 0, 1.2, 100, -5.0, 5.0);
                    }
                }
            }
        }
    }

    void AnalysisTiming::run(DataBatch& batch) {
        if (batch.tracks.empty()) return;
        auto& geo = GeometryManager::getInstance();

        size_t n_tracks = batch.tracks.size();
        size_t n_threads = thread_pool->getThreadCount();
        size_t chunk = (n_tracks + n_threads - 1) / n_threads;

        std::vector<std::future<void>> futures;

        for (size_t i = 0; i < n_tracks; i += chunk) {
            size_t end = std::min(i + chunk, n_tracks);
            
            futures.push_back(thread_pool->submit([&, i, end]() {
                size_t loc_valid = 0;
                size_t loc_skipped = 0;
                auto& pm = PlotManager::getInstance();
                // Only touched when time_walk_correction_ is set - see its
                // docs. Accumulated locally per thread-chunk, merged into
                // the shared members once at the end of this lambda (NOT
                // per-cluster/per-pair - #deferred_mtx_ would be
                // hammered otherwise).
                std::map<std::string, TimeWalkFitAccum> loc_fit;
                std::vector<PairRecord> loc_pairs;

                for (size_t t_idx = i; t_idx < end; ++t_idx) {
                    auto const& track = batch.tracks[t_idx];

                    // 1. Check required planes (reference planes like RD53B) -
                    // matches corry's own foundRequiredAssocCluster():
                    // required_detectors_ is a separate, explicit list from
                    // det_names_ (the timing detectors themselves).
                    bool req_met = true;
                    for (auto const& req : required_detectors_) {
                        if (!track->hasDetector(req)) { req_met = false; break; }
                    }
                    if (!req_met) continue;

                    // 2. Bucket this track's associated clusters by detector
                    // name in one pass. A vector per detector rather than
                    // a single cluster, to faithfully support a track
                    // associated to more than one cluster on the same
                    // detector - matches corry's own cluster_vec_i loop,
                    // which doesn't break after the first match either (only
                    // its separate, earlier validity-gating loop does).
                    std::map<std::string, std::vector<std::shared_ptr<Cluster>>> by_det;
                    for (auto const& c : track->getAssociatedClusters()) {
                        by_det[c->detectorID()].push_back(c);
                    }

                    // 3. corry requires EVERY configured detector to have at
                    // least one associated cluster before a track
                    // contributes anything at all ("Check that the track has
                    // a hit in all DUTs" in its run()) - matched here
                    // exactly, not relaxed to a per-pair requirement.
                    bool valid_track = true;
                    for (auto const& name : det_names_) {
                        if (by_det.find(name) == by_det.end()) { valid_track = false; break; }
                    }
                    if (!valid_track) continue;

                    loc_valid++;

                    auto raw_timestamp = [](const std::shared_ptr<Cluster>& c) { return c->earliestPixelTimestamp(); };

                    // Optional, Warlock-only, off by default: t_corr = t_raw
                    // - (a/charge + b) or t_raw - (a/sqrt(charge) + b).
                    // corry's own module never applies a time-walk
                    // correction at all (its run() has a literal
                    // "// MISSING --- Time walk correction" comment) - this
                    // is a no-op matching corry exactly unless configured.
                    auto corrected_timestamp = [this](const std::string& det, double t_raw, double charge) {
                        if (charge <= 0 || time_walk_method_ == "none") return t_raw;
                        double a = time_walk_a_.at(det);
                        double b = time_walk_b_.at(det);
                        if (time_walk_method_ == "inverse_sqrt_charge") return t_raw - (a / std::sqrt(charge) + b);
                        if (time_walk_method_ == "inverse_charge") return t_raw - (a / charge + b);
                        return t_raw;
                    };

                    // 4. For every detector i and every one of its
                    // associated clusters (usually exactly one), fill
                    // charge_toa_i against the RAW (uncorrected) earliest-
                    // pixel timestamp - matches corry exactly, which never
                    // applies time-walk correction here either - then pair
                    // against every detector j>i's associated cluster(s).
                    // Mirrors corry's nested cluster_vec_i/cluster_vec_j
                    // loop structure, including the rare multi-cluster-per-
                    // detector case.
                    for (size_t d1_idx = 0; d1_idx < det_names_.size(); ++d1_idx) {
                        const std::string& name1 = det_names_[d1_idx];
                        auto& det1 = geo.getDetector(name1);
                        auto const& clusters1 = by_det.at(name1);
                        int n_pairs = 0;
                        // Per-die (TOP/BOTTOM) breakdown of ncat_/npairs_ -
                        // "" would double-count the combined total already
                        // filled below, so it's deliberately skipped when
                        // flushing these maps.
                        std::map<std::string, int> ncat_by_die, npairs_by_die;

                        for (auto const& c1 : clusters1) {
                            std::string die1 = det1.dieOf(c1->column(), c1->row());
                            std::string label1 = die1.empty() ? name1 : name1 + "_" + die1;
                            ncat_by_die[die1]++;

                            double a1 = c1->charge();
                            double t1_raw = raw_timestamp(c1);
                            double t1 = corrected_timestamp(name1, t1_raw, a1);
                            pm.fill2D(getName(), "charge_toa_" + name1, a1, t1_raw);
                            if (!die1.empty()) pm.fill2D(getName(), "charge_toa_" + label1, a1, t1_raw);

                            int n_pairs_c1 = 0;
                            for (size_t d2_idx = d1_idx + 1; d2_idx < det_names_.size(); ++d2_idx) {
                                const std::string& name2 = det_names_[d2_idx];
                                auto& det2 = geo.getDetector(name2);
                                for (auto const& c2 : by_det.at(name2)) {
                                    // Skip cluster pairs sampled by different
                                    // CAEN trigger groups - not physically
                                    // comparable timing (see
                                    // require_same_trigger_group_'s docs).
                                    // triggerGroupOf() returns -1 for a
                                    // channel outside 0-17 (e.g. a non-CAEN
                                    // detector's default -1) - explicitly
                                    // excluded from matching (not just
                                    // compared with !=), so two "unknown"
                                    // clusters never spuriously pass as
                                    // same-group.
                                    if (require_same_trigger_group_) {
                                        int g1 = triggerGroupOf(c1->seedChannel());
                                        int g2 = triggerGroupOf(c2->seedChannel());
                                        if (g1 < 0 || g2 < 0 || g1 != g2) {
                                            ++loc_skipped;
                                            continue;
                                        }
                                    }
                                    std::string die2 = det2.dieOf(c2->column(), c2->row());
                                    std::string label2 = die2.empty() ? name2 : name2 + "_" + die2;
                                    double a2 = c2->charge();
                                    double t2_raw = raw_timestamp(c2);

                                    if (time_walk_correction_) {
                                        // The correction isn't known until
                                        // finalize() has fit
                                        // time_walk_fit_accum_ from every
                                        // batch - stash the raw inputs and
                                        // fill dtoa_*/timewalk_* there
                                        // instead of here.
                                        loc_pairs.push_back({name1, die1, name2, die2, a1, t1_raw, a2, t2_raw});

                                        // Fit input for BOTH detectors in
                                        // this pair - local_delta is each
                                        // one's own raw ToA minus the OTHER
                                        // side's, so both are small/well-
                                        // scaled regardless of how far into
                                        // the run this event happened (see
                                        // #time_walk_correction_'s docs).
                                        // Every pair a detector appears in
                                        // (as either name1 or name2)
                                        // contributes - not just pairs where
                                        // it happens to be name1. Guarded by
                                        // a real floor, not just `> 0`: a
                                        // near-zero (but positive) charge
                                        // outlier makes x=1/charge blow up
                                        // to O(1e11)+, and that single point
                                        // then dominates the OLS sums via
                                        // extreme leverage.
                                        constexpr double kMinFitCharge = 1e-3;
                                        if (a1 > kMinFitCharge) {
                                            double x1 = (time_walk_method_ == "inverse_sqrt_charge") ? 1.0 / std::sqrt(a1) : 1.0 / a1;
                                            double local_delta1 = t1_raw - t2_raw;
                                            auto& fa1 = loc_fit[name1];
                                            fa1.sum_x += x1; fa1.sum_x2 += x1 * x1;
                                            fa1.sum_t += local_delta1; fa1.sum_xt += x1 * local_delta1;
                                            fa1.n += 1;
                                        }
                                        if (a2 > kMinFitCharge) {
                                            double x2 = (time_walk_method_ == "inverse_sqrt_charge") ? 1.0 / std::sqrt(a2) : 1.0 / a2;
                                            double local_delta2 = t2_raw - t1_raw;
                                            auto& fa2 = loc_fit[name2];
                                            fa2.sum_x += x2; fa2.sum_x2 += x2 * x2;
                                            fa2.sum_t += local_delta2; fa2.sum_xt += x2 * local_delta2;
                                            fa2.n += 1;
                                        }
                                    } else {
                                        double t2 = corrected_timestamp(name2, t2_raw, a2);
                                        double delta_t = t1 - t2;
                                        std::string pair = name1 + "_vs_" + name2;

                                        pm.fill1D(getName(), "dtoa_" + pair, delta_t);
                                        pm.fill2D(getName(), "timewalk_" + name1 + "_in_" + pair, a1, delta_t);

                                        // Die-qualified variant, only if this
                                        // pair involves at least one split
                                        // detector (label1/label2 differs from
                                        // name1/name2) - matches c1's/c2's own
                                        // die independently, so a TOP-vs-BOTTOM
                                        // cross-die pairing is tracked distinctly
                                        // from TOP-vs-TOP or combined.
                                        if (!die1.empty() || !die2.empty()) {
                                            std::string pair_variant = label1 + "_vs_" + label2;
                                            pm.fill1D(getName(), "dtoa_" + pair_variant, delta_t);
                                            pm.fill2D(getName(), "timewalk_" + label1 + "_in_" + pair_variant, a1, delta_t);
                                        }
                                    }

                                    ++n_pairs;
                                    ++n_pairs_c1;
                                }
                            }
                            npairs_by_die[die1] += n_pairs_c1;
                        }
                        pm.fill1D(getName(), "ncat_" + name1, static_cast<double>(clusters1.size()));
                        pm.fill1D(getName(), "npairs_" + name1, static_cast<double>(n_pairs));
                        for (auto const& [die1, cnt] : ncat_by_die) {
                            if (die1.empty()) continue;
                            pm.fill1D(getName(), "ncat_" + name1 + "_" + die1, static_cast<double>(cnt));
                        }
                        for (auto const& [die1, cnt] : npairs_by_die) {
                            if (die1.empty()) continue;
                            pm.fill1D(getName(), "npairs_" + name1 + "_" + die1, static_cast<double>(cnt));
                        }
                    }
                }
                valid_timing_tracks_ += loc_valid;
                pairs_skipped_group_mismatch_ += loc_skipped;

                if (time_walk_correction_) {
                    std::lock_guard<std::mutex> lock(deferred_mtx_);
                    for (auto const& [det, fa] : loc_fit) {
                        auto& shared = time_walk_fit_accum_[det];
                        shared.sum_x += fa.sum_x; shared.sum_x2 += fa.sum_x2;
                        shared.sum_t += fa.sum_t; shared.sum_xt += fa.sum_xt;
                        shared.n += fa.n;
                    }
                    pair_records_.insert(pair_records_.end(),
                                          std::make_move_iterator(loc_pairs.begin()),
                                          std::make_move_iterator(loc_pairs.end()));
                }
            }));
        }
        for (auto& f : futures) f.get();
        total_tracks_ += n_tracks;
    }

    void AnalysisTiming::finalize() {
        WR_LOG(STATUS, "Timing Analysis Finished.");
        WR_LOG(INFO, "Tracks with a cluster on every configured detector: " + std::to_string(valid_timing_tracks_) +
                     " / " + std::to_string(total_tracks_));

        if (time_walk_method_ != "none") {
            WR_LOG(INFO, "Time-walk correction applied: method=" + time_walk_method_);
        }

        if (time_walk_correction_) {
            auto& pm = PlotManager::getInstance();

            // Closed-form ordinary least squares per detector: minimizes
            // sum((local_delta - a*x - b)^2) over that detector's pooled
            // (x, local_delta) points (see #time_walk_correction_ for why
            // local_delta, not raw t_raw). `b` here is the fitted
            // intercept - logged for visibility (it's the average
            // propagation-delay/latency offset vs this detector's various
            // partners) but deliberately NOT used in the correction itself;
            // see fittedCorrect() below.
            std::map<std::string, std::pair<double, double>> fitted_a_meanx;
            for (auto const& [det, fa] : time_walk_fit_accum_) {
                double n = static_cast<double>(fa.n);
                double denom = n * fa.sum_x2 - fa.sum_x * fa.sum_x;
                if (fa.n < 2 || std::abs(denom) < 1e-30) {
                    WR_LOG(WARNING, det + ": not enough spread in charge to fit a time-walk correction (" +
                                         std::to_string(fa.n) + " point(s)) - left uncorrected for this run.");
                    fitted_a_meanx[det] = {0.0, 0.0};
                    continue;
                }
                double a = (n * fa.sum_xt - fa.sum_x * fa.sum_t) / denom;
                double b = (fa.sum_t - a * fa.sum_x) / n;
                double mean_x = fa.sum_x / n;
                fitted_a_meanx[det] = {a, mean_x};
                WR_LOG(STATUS, "AnalysisTiming: fitted time-walk for " + det + " (method=" + time_walk_method_ +
                                    ", n=" + std::to_string(fa.n) + "): a=" + std::to_string(a) +
                                    ", intercept=" + std::to_string(b) + " (not applied - see docs), mean(x)=" +
                                    std::to_string(mean_x));
            }

            // t_corr = t_raw - a*(x - mean_x): zero correction at this
            // detector's own mean charge, removing only the charge-
            // dependent slope - see #time_walk_correction_'s docs for why
            // the fit's own intercept is deliberately not subtracted here.
            auto fittedCorrect = [&](const std::string& det, double t_raw, double charge) {
                // Same floor as the fit's own input filter above - a
                // near-zero-charge cluster wasn't part of the fit and
                // shouldn't get an extreme x=1/charge correction applied
                // to it either.
                constexpr double kMinFitCharge = 1e-3;
                if (charge <= kMinFitCharge) return t_raw;
                auto it = fitted_a_meanx.find(det);
                if (it == fitted_a_meanx.end()) return t_raw;
                double x = (time_walk_method_ == "inverse_sqrt_charge") ? 1.0 / std::sqrt(charge) : 1.0 / charge;
                return t_raw - it->second.first * (x - it->second.second);
            };

            // Replay every pair run() deferred (see #pair_records_) now
            // that the fit is known - mirrors run()'s own dtoa_*/timewalk_*
            // fill logic exactly, just fed from stored records instead of
            // live Cluster objects, and using the fitted a/b instead of
            // config-supplied ones.
            for (auto const& rec : pair_records_) {
                double t1 = fittedCorrect(rec.name1, rec.t1_raw, rec.charge1);
                double t2 = fittedCorrect(rec.name2, rec.t2_raw, rec.charge2);
                double delta_t = t1 - t2;
                std::string pair = rec.name1 + "_vs_" + rec.name2;

                pm.fill1D(getName(), "dtoa_" + pair, delta_t);
                pm.fill2D(getName(), "timewalk_" + rec.name1 + "_in_" + pair, rec.charge1, delta_t);

                if (!rec.die1.empty() || !rec.die2.empty()) {
                    std::string label1 = dieLabel(rec.name1, rec.die1);
                    std::string label2 = dieLabel(rec.name2, rec.die2);
                    std::string pair_variant = label1 + "_vs_" + label2;
                    pm.fill1D(getName(), "dtoa_" + pair_variant, delta_t);
                    pm.fill2D(getName(), "timewalk_" + label1 + "_in_" + pair_variant, rec.charge1, delta_t);
                }
            }
            WR_LOG(STATUS, "AnalysisTiming: filled dtoa_*/timewalk_* plots for " + std::to_string(pair_records_.size()) +
                                " pair(s) using the fitted time-walk correction.");
        }

        if (require_same_trigger_group_) {
            WR_LOG(INFO, "Cluster pairs skipped for trigger-group mismatch: " + std::to_string(pairs_skipped_group_mismatch_));
            if (pairs_skipped_group_mismatch_ > 0 && valid_timing_tracks_ > 0 &&
                pairs_skipped_group_mismatch_ >= valid_timing_tracks_ * (det_names_.size() * (det_names_.size() - 1) / 2)) {
                WR_LOG(WARNING, "Every cluster pair was skipped for trigger-group mismatch - if this is unexpected, "
                                 "check that the input file's clusters carry real channel info (re-saved after "
                                 "channel tracking was added to Tracking4D's cluster schema), or set "
                                 "require_same_trigger_group=false to disable this check.");
            }
        }
    }

    REGISTER_MODULE(AnalysisTiming)
}