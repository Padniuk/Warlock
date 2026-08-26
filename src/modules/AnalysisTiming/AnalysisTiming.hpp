/**
 * @file AnalysisTiming.hpp
 * @brief Cross-detector time-of-arrival correlation and optional time-walk correction.
 */

#ifndef WARLOCK_ANALYSISTIMING_HPP
#define WARLOCK_ANALYSISTIMING_HPP

#include "core/Module.hpp"
#include <atomic>
#include <string>
#include <vector>
#include <map>
#include <mutex>

namespace framework {
    /**
     * @brief Mirrors corryvreckan's AnalysisTiming module: for every track
     * carrying an associated cluster on ALL of #det_names_ (and any
     * #required_detectors_), fills each detector's own charge-vs-ToA map
     * plus a pairwise delta-ToA histogram for every detector pair, using
     * the earliest pixel timestamp per cluster. Time-walk correction
     * (#time_walk_method_) is a Warlock-only addition, off by default -
     * corry's own module never applies one (its run() has a literal
     * "MISSING --- Time walk correction" comment).
     */
    class AnalysisTiming : public Module {
    public:
        AnalysisTiming(Configuration cfg, Configuration g_cfg, ThreadPool* pool);

        void initialize() override;
        void run(DataBatch& batch) override;
        void finalize() override;

    private:
        std::vector<std::string> det_names_;         ///< Detectors correlated against each other; a track needs a cluster on every one to count.
        std::vector<std::string> required_detectors_; ///< Additional detectors (e.g. a reference plane) a track must also have a cluster on, but that aren't themselves correlated.

        /// Time-walk correction: t_corr = t_raw - (a/charge + b) or
        /// t_raw - (a/sqrt(charge) + b), matching corry's parametrization.
        std::string time_walk_method_; ///< "none" (default), "inverse_charge", or "inverse_sqrt_charge".
        std::map<std::string, double> time_walk_a_; ///< Per-detector `a` coefficient, keyed by detector name. Ignored (overwritten at finalize()) when #time_walk_correction_ is set.
        std::map<std::string, double> time_walk_b_; ///< Per-detector `b` coefficient, keyed by detector name. Ignored (overwritten at finalize()) when #time_walk_correction_ is set.

        /// Warlock-only addition: if set, `a`/`b` are NOT read from config at
        /// all - instead fit per detector from this run's own data, then
        /// applied to that detector's dtoa_*/timewalk_* histograms.
        ///
        /// The fit target is deliberately NOT `t_raw` itself: a cluster's
        /// raw `earliestPixelTimestamp()` is an absolute value spanning the
        /// whole run, so pooling (charge, t_raw) across every event mostly
        /// measures when in the run the trigger happened rather than any
        /// charge dependence, swamping the genuinely small (ns-scale)
        /// charge-walk signal. Instead, for every cluster pair that passes
        /// run()'s cuts, each detector D in the pair contributes one point
        /// `(f(charge_D), t_D_raw - t_other_raw)` to its own fit - the
        /// pair's own delta-ToA, already small and well-scaled since both
        /// timestamps come from the same physical event. A detector
        /// contributes from every pair it appears in (as either side).
        ///
        /// Correction applied: `t_corr = t_raw - a*(x - mean_x)`, i.e.
        /// anchored so a cluster at this detector's own mean charge gets no
        /// correction - only deviations from typical charge get shifted.
        /// The fit's own intercept is intentionally not subtracted: it
        /// reflects genuine propagation-delay/trigger-latency differences
        /// between this detector and its partners, which is real physics
        /// to preserve, not a time-walk artifact to remove.
        ///
        /// Requires #time_walk_method_ to be set to a real method (not
        /// "none") - if it's still "none" when this is set, a WARNING is
        /// logged in the constructor and this flag is forced back to false.
        /// Since the fit isn't known until every batch has been seen,
        /// dtoa_*/timewalk_* histograms are filled at finalize() instead of
        /// run() in this mode - #pair_records_ keeps the lean per-pair data
        /// needed for that replay.
        bool time_walk_correction_{false};

        /// Ordinary-least-squares running sums for one detector's own
        /// `local_delta = a*x + b` fit (x = f(charge), local_delta = this
        /// detector's own raw ToA minus its pair partner's, see
        /// #time_walk_correction_ for why NOT plain t_raw) - O(1) memory
        /// regardless of cluster count, standard closed-form linear
        /// regression (only these five sums are ever needed). `sum_t`/
        /// `sum_xt` accumulate local_delta, not t_raw, despite the name
        /// (kept short/generic since the struct is a plain x/y OLS
        /// accumulator either way).
        struct TimeWalkFitAccum {
            double sum_x{0.0}, sum_x2{0.0}, sum_t{0.0}, sum_xt{0.0};
            size_t n{0};
        };
        std::map<std::string, TimeWalkFitAccum> time_walk_fit_accum_; ///< Keyed by detector name (#det_names_ entries) - guarded by #deferred_mtx_.

        /// One entry per (cluster1, cluster2) pair that passed every run()
        /// cut (required detectors, trigger-group match) - only populated
        /// when #time_walk_correction_ is set, since the correction (and
        /// hence the dtoa/timewalk value to fill) isn't known until
        /// finalize() has fit #time_walk_fit_accum_. Kept lean (a handful
        /// of small strings and doubles per pair) to stay cheap at the
        /// track counts this module runs over.
        struct PairRecord {
            std::string name1, die1, name2, die2; ///< die1/die2 are "" (unsplit) or one of the detector's own DetectorGeo::dieLabels() (e.g. "TOP"/"BOTTOM").
            double charge1{0.0}, t1_raw{0.0}, charge2{0.0}, t2_raw{0.0};
        };
        std::vector<PairRecord> pair_records_; ///< Guarded by #deferred_mtx_.
        std::mutex deferred_mtx_; ///< Guards #time_walk_fit_accum_/#pair_records_ - locked once per run() thread-chunk (merging that chunk's local accumulation), never per-pair.

        std::atomic<size_t> total_tracks_{0};         ///< Every track seen in run(), across all batches.
        std::atomic<size_t> valid_timing_tracks_{0};  ///< Subset of #total_tracks_ with a cluster on every configured detector.

        /// If set (default true), a cluster pair is only compared (all
        /// charge_toa/dtoa/timewalk/npairs fills) when both clusters'
        /// Cluster::seedChannel()-derived trigger groups match - comparing
        /// timing across two CAEN channels sampled by different trigger
        /// groups on the digitizer isn't physically meaningful (each group
        /// free-runs independently). Warlock-only, corry has no concept of
        /// this at all (see Cluster.hpp's triggerGroupOf()). Applied
        /// pairwise but correctly generalizes to N detectors via
        /// transitivity - if A-B share group X and A-C share group X, B-C
        /// automatically share group X too.
        bool require_same_trigger_group_;
        std::atomic<size_t> pairs_skipped_group_mismatch_{0}; ///< Diagnostic counter, reported in finalize().
    };
}
#endif
