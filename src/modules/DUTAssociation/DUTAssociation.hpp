#ifndef WARLOCK_DUTASSOCIATION_HPP
#define WARLOCK_DUTASSOCIATION_HPP

#include "core/Module.hpp"
#include <string>
#include <vector>
#include <map>
#include <mutex>
#include <atomic>

namespace framework {
    /**
     * @brief Matches fitted tracks to clusters on one or more DUTs, mirroring
     * corryvreckan's DUTAssociation module.
     *
     * For every track and every configured target detector, associates
     * every cluster within an elliptical spatial cut (@ref spatial_cut_x_ /
     * @ref spatial_cut_y_) and an absolute time cut (@ref time_cut_abs_) -
     * not just the closest one, matching corry's own behavior. Bound to
     * either a single detector (`name`) or every detector of a given
     * `type`, one instance per config block. For a two-die CAEN board
     * (DetectorGeo::hasDieSplit()), every plot is filled both combined and
     * split into `/TOP` and `/BOTTOM` subpaths.
     *
     * Candidate clusters are restricted to the die the TRACK itself
     * predicts (DetectorGeo::getColumn()/getRow() of its own local
     * intercept, then dieOf()) - not every cluster on the detector
     * regardless of die (Warlock-only fix, 2026-09-01). Without this, a
     * track near the physical gap between two dies gets tested against
     * both dies' boundary-row clusters, each an independent chance to
     * fall inside the ellipse cut - visible as efficiency/charge maps
     * "bleeding" further into the gap between two dies than past the true
     * outer edge of the whole assembly, where only one die's clusters are
     * ever nearby. No-op for a detector with no die split at all.
     */
    class DUTAssociation : public Module {
    public:
        DUTAssociation(Configuration cfg, Configuration g_cfg, ThreadPool* pool);

        void initialize() override;
        void run(DataBatch& batch) override;
        void finalize() override;

    private:
        // Which per-die plots (if any) this detector splits into - decided
        // via the single shared DetectorGeo::dieLabels() (GeometryManager.hpp),
        // not re-derived here, so AnalysisEfficiency/AnalysisTiming's
        // identical decision can never drift from this one. Read once per
        // detector from initialize() and cached in #die_labels_by_det_ -
        // NOT recomputed per-track from the multi-threaded run() loop.
        std::map<std::string, std::vector<std::string>> die_labels_by_det_; ///< Cached DetectorGeo::dieLabels() result per target detector, populated in initialize(). Empty entry means "doesn't split".

        std::string dut_name_;    // explicit single-detector binding ("name")
        std::string type_filter_; // multi-detector binding ("type")
        std::vector<std::string> target_detectors_;
        double spatial_cut_x_;
        double spatial_cut_y_;
        double time_cut_abs_;
        bool use_cluster_centre_;

        std::atomic<size_t> total_tracks_{0};
        std::atomic<size_t> total_associations_{0};
        std::atomic<size_t> tracks_with_assoc_{0};
        std::atomic<size_t> cut_spatial_{0};
        std::atomic<size_t> cut_timing_{0};

        // Count of every (track, cluster) pair examined (associated or cut)
        // - corry's "num_cluster", used to normalize hCutHisto into a
        // fraction at finalize() the same way corry does
        // (hCutHisto->Scale(1 / double(num_cluster))). Keyed by plot
        // subpath, not just detector name: a two-die board
        // (DetectorGeo::hasDieSplit()) gets three independent entries -
        // "<det_name>" (both dies combined, unchanged from before the
        // split), "<det_name>/TOP", "<det_name>/BOTTOM" - each normalizing
        // its own hCutHisto against its own examined-pair count.
        std::map<std::string, size_t> num_examined_by_det_;
        std::mutex counter_mtx_;
    };
}
#endif
