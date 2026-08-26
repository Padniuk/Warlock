/**
 * @file Prealignment.hpp
 * @brief Coarse-shifts each detector's position to align it with the
 * reference plane, using the peak of its cluster correlation against the
 * reference.
 */

#ifndef WARLOCK_PREALIGNMENT_HPP
#define WARLOCK_PREALIGNMENT_HPP

#include "core/Module.hpp"
#include <string>
#include <vector>

namespace framework {
    /**
     * @brief Correlates every detector's clusters against the reference
     * plane's, then shifts each detector's position by the correlation
     * peak/mean (per `method`) - a single-pass, non-iterative first-order
     * alignment, typically run before AlignmentTrackChi2 or
     * AlignmentMillepede. Mirrors corryvreckan's Prealignment module,
     * including instantiating (and histogramming) the reference detector
     * itself, matching it against its own clusters - only the actual
     * position shift is skipped for it.
     */
    class Prealignment : public Module {
    public:
        Prealignment(Configuration cfg, Configuration g_cfg, ThreadPool* pool);

        void initialize() override;
        void run(DataBatch& batch) override;
        /// Computes and applies each non-reference, non-fixed detector's
        /// shift from its accumulated correlation histograms, then saves
        /// the updated geometry to `detectors_file_updated`.
        void finalize() override;

    private:
        std::string method_;       ///< `"mean"`, `"maximum"`, or `"gauss_fit"` - which statistic of the correlation histogram to shift by.
        std::string type_filter_;  ///< Only detectors of this `type` are aligned; empty = all pixel detectors.

        double range_abs_;           ///< Correlation histogram half-range, in mm.
        double damping_factor_;      ///< Scales the computed shift before applying it (1.0 = full shift).
        double max_correlation_rms_; ///< RMS threshold above which a warning is logged (shift is still applied).
        double time_cut_abs_;        ///< Cluster-pair time cut, in ns.
        std::vector<std::string> fixed_planes_; ///< Detector names excluded from shifting (histogrammed like any other).
    };
}
#endif
