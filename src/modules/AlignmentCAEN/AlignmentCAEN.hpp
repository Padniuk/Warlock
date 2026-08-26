/**
 * @file AlignmentCAEN.hpp
 * @brief XY translation and Z-rotation alignment for CAEN DUT boards, driven purely by DUTAssociation-matched tracks.
 */

#ifndef WARLOCK_ALIGNMENTCAEN_HPP
#define WARLOCK_ALIGNMENTCAEN_HPP

#include "core/Module.hpp"
#include <Eigen/Dense>
#include <map>
#include <mutex>
#include <string>
#include <vector>

namespace framework {
    /**
     * @brief Dedicated alignment module driven purely by DUTAssociation's
     * output (which tracks have an associated cluster on this detector)
     * plus each track's own fitted intercept. No efficiency map,
     * charge/cluster-size/ToA profiling, or fake-rate accounting is
     * computed here - that separation of "align" (this module,
     * PrealignmentCAEN, AlignmentTrackChi2, AlignmentMillepede) from
     * "measure" (AnalysisEfficiency, AnalysisTiming) is deliberate.
     *
     * Config resolution and track cut-flow match AnalysisEfficiency's own
     * (same "name"/"type" binding, same chi2/DUT-edge/masked-pixel/
     * required-associated-cluster gates), so a config block's alignment-
     * relevant keys can move from [AnalysisEfficiency] to [AlignmentCAEN]
     * without changing which tracks count toward the fit.
     */
    class AlignmentCAEN : public Module {
    public:
        AlignmentCAEN(Configuration cfg, Configuration g_cfg, ThreadPool* pool);

        void initialize() override;
        void run(DataBatch& batch) override;
        void finalize() override;

    private:
        std::string dut_name_;    ///< Explicit single-detector binding ("name").
        std::string type_filter_; ///< Multi-detector binding ("type").
        std::vector<std::string> target_detectors_;

        double chi2_cut_;                            ///< Track chi2/ndof cut, matching AnalysisEfficiency's `chi2ndof_cut`.
        std::vector<std::string> required_detectors_; ///< `require_associated_cluster_on` - a track must have an associated cluster on every one of these too, or it's dropped before it can count toward either fit.
        int masked_pixel_distance_cut_;
        double spatial_cut_sensoredge_;

        bool update_xy_translation_;
        bool update_z_orientation_;
        int z_rotation_iterations_;
        /// Manual override on top of the automatic per-pixel check below:
        /// 0 (default) uses every pixel the per-pixel check accepts;
        /// positive/negative restricts to the top/bottom row outright,
        /// for a case an operator already knows is bad even if its
        /// statistics happen to look locally plausible.
        int z_rotation_row_;
        double min_quadrant_entries_; ///< Minimum weighted matched-intercept count a pixel needs before its centroid is trusted at all.
        /// Maximum allowed distance (as a fraction of the smaller pitch axis)
        /// between a pixel's measured intercept centroid and its nominal
        /// position before the pixel is excluded from both the translation
        /// pool and the rotation fit. Applied per pixel rather than per row,
        /// so single-row boards get the same dead-pixel protection as
        /// multi-row ones.
        double max_pixel_centroid_deviation_;
        double bin_size_x_, bin_size_y_; ///< Local-frame bin size (mm) for the matched-intercept grid pixelCentroid() reads from.

        struct TrackRecord {
            Eigen::Vector3d global_pred;
            bool matched;
        };
        /// Per-detector accumulator - one per entry in #target_detectors_.
        struct DetAccum {
            int n_bins_x{0}, n_bins_y{0};
            double size_x{0}, size_y{0};
            /// Matched-track-intercept counts, local frame - rebuilt from
            /// #rotation_records against the CURRENT geometry each time
            /// finalize() needs it (see rebinLocal() in the .cpp), since
            /// a rotation/translation correction changes which bin (and
            /// which pixel) an intercept falls into.
            Eigen::MatrixXd eff_num_local;
            /// One entry per track that survived the full cut-flow
            /// (DUT-bounds/chi2/isolation/masked-pixel/required-cluster),
            /// whether or not it ended up matched - global-frame, so
            /// finalize() can re-derive local_pred against an updated
            /// R_inv/position on every rotation-fit iteration instead of
            /// only ever seeing this run's original geometry.
            std::vector<TrackRecord> rotation_records;
        };
        std::map<std::string, DetAccum> accum_;
        std::mutex acc_mtx_; ///< Guards merging each thread chunk's local records into #accum_.

        /// Applies the same track selection cut-flow AnalysisEfficiency
        /// uses (DUT bounds, chi2, isolation, masked-pixel proximity,
        /// required associated clusters) and appends one TrackRecord per
        /// accepted track - no plot filling, no efficiency/profile
        /// accumulation, just the raw (position, matched) pairs the
        /// finalize()-time fits need.
        void runTrackLoop(DataBatch& batch, const std::string& det_name);
    };
}
#endif
