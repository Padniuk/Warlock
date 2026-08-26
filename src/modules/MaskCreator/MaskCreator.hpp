/**
 * @file MaskCreator.hpp
 * @brief Identifies and masks noisy/dead pixels from per-detector occupancy.
 */

#ifndef WARLOCK_MASKCREATOR_HPP
#define WARLOCK_MASKCREATOR_HPP

#include "core/Module.hpp"
#include <Eigen/Dense>
#include <atomic>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <vector>

namespace framework {
    /**
     * @brief Builds a per-detector pixel occupancy map across the run, then
     * flags outlier pixels as masked using one of two methods (see
     * `method` config key): `"frequency"` (a pixel is masked if its hit
     * count exceeds `frequency_cut` times the detector's mean occupancy)
     * or `"localdensity"` (a kernel-density comparison against each
     * pixel's own neighborhood, matching corry's `Ignore` module). Mirrors
     * corryvreckan's MaskCreator, built on Warlock's batch-parallel
     * architecture instead.
     */
    class MaskCreator : public Module {
    public:
        /**
         * @param cfg Module-specific configuration.
         * @param g_cfg Global framework configuration.
         * @param pool The thread pool shared by the whole pipeline.
         */
        MaskCreator(Configuration cfg, Configuration g_cfg, ThreadPool* pool);

        void initialize() override;
        void run(DataBatch& batch) override;
        void finalize() override;

    private:
        /// Computes the `"localdensity"` mask for one detector: for each
        /// occupied pixel, compares its hit count against the mean of its
        /// own `density_bandwidth`-pixel neighborhood, flagging it if the
        /// local significance exceeds `sigma_above_avg_max` or its
        /// per-event rate exceeds `rate_max`. Sharded by row across the
        /// thread pool. Also fills the `density`/`local_significance`
        /// diagnostic plots from the same computation.
        void calculateDensityMask(const std::string& det_name, const Eigen::MatrixXd& matrix,
                                   std::vector<std::pair<int, int>>& masked_indices, Eigen::MatrixXd& mask_display);

        std::string type_filter_;   ///< Only detectors of this `type` are masked; empty = all pixel detectors.
        std::string method_;        ///< `"frequency"` or `"localdensity"`.
        double frequency_cut_;      ///< `"frequency"` method threshold, as a multiple of the mean.
        int density_bandwidth_;     ///< `"localdensity"` method neighborhood half-size, in pixels.
        double sigma_max_;          ///< `"localdensity"` method significance threshold.
        double rate_max_;           ///< `"localdensity"` method per-event rate threshold.
        bool mask_dead_pixels_;     ///< Whether never-hit pixels are also masked.

        std::atomic<size_t> total_events_{0};

        std::map<std::string, Eigen::MatrixXd> occupancy_matrices_;       ///< Per-detector hit counts, accumulated across every batch.
        std::map<std::string, std::unique_ptr<std::mutex>> detector_mutexes_; ///< One lock per detector, guarding its own #occupancy_matrices_ entry.

        /// Pixels already masked via a pre-existing mask_file, preserved and
        /// merged with newly-found ones - matches corry pre-filling its mask
        /// map from previously known masks instead of overwriting them.
        std::map<std::string, std::set<std::pair<int, int>>> existing_masks_;
    };
}

#endif // WARLOCK_MASKCREATOR_HPP
