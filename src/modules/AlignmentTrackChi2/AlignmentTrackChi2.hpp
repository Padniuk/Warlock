/**
 * @file AlignmentTrackChi2.hpp
 * @brief Iterative track-based alignment by minimizing per-track fit chi2.
 */

#ifndef WARLOCK_ALIGNMENTTRACKCHI2_HPP
#define WARLOCK_ALIGNMENTTRACKCHI2_HPP

#include "core/Module.hpp"
#include <Eigen/Dense>
#include <mutex>
#include <vector>
#include <string>

namespace framework {
    /**
     * @brief For each non-reference, non-fixed detector in turn, runs an
     * Eigen Levenberg-Marquardt fit (via AlignmentFunctor, defined in the
     * .cpp) that perturbs the detector's position/orientation and refits
     * every track touching it, minimizing sqrt(chi2) per track. Repeats for
     * #iterations_ passes over all detectors. Architecturally similar to
     * corryvreckan's own AlignmentTrackChi2 (accumulate tracks, then iterate
     * per-detector minimization), but uses a different numerical optimizer,
     * so convergence is not bit-exact - agreement is judged by final
     * detector position, not by intermediate values.
     */
    class AlignmentTrackChi2 : public Module {
    public:
        AlignmentTrackChi2(Configuration cfg, Configuration g_cfg, ThreadPool* pool);

        void initialize() override;
        void run(DataBatch& batch) override;
        void finalize() override;

        double getProgress() const override {
            if (max_tracks_ > 0) return static_cast<double>(total_evaluated_tracks_) / max_tracks_;
            return -1.0;
        }

    private:
        int iterations_;             ///< Number of full passes over all detectors.
        double max_chi2ndof_;        ///< Chi2/ndof threshold used by #prune_tracks_.
        bool prune_tracks_;          ///< If set, tracks above #max_chi2ndof_ are dropped before alignment.
        bool align_position_;        ///< Whether X/Y position are free alignment parameters.
        bool align_orientation_;     ///< Whether X/Y/Z rotation are free alignment parameters.
        size_t max_tracks_;          ///< Framework-wide track target (corry's `number_of_tracks`), used only for #getProgress().
        std::string type_filter_;    ///< If non-empty, only detectors of this type are aligned.
        std::vector<std::string> fixed_planes_; ///< Detectors excluded from alignment despite being non-reference.

        std::vector<std::shared_ptr<Track>> global_tracks_; ///< Tracks accumulated across all batches, aligned against in finalize().
        size_t total_evaluated_tracks_{0}; ///< Count of every track seen in run(), including pruned ones.
        std::mutex track_mtx_; ///< Guards #global_tracks_ and #total_evaluated_tracks_ across concurrent run() calls.
    };
}
#endif
