/**
 * @file AnalysisEfficiency.hpp
 * @brief DUT tracking efficiency, charge/cluster profile maps, and fake-hit rate.
 */

#ifndef WARLOCK_ANALYSISEFFICIENCY_HPP
#define WARLOCK_ANALYSISEFFICIENCY_HPP

#include "core/Module.hpp"
#include <Eigen/Dense>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace framework {
    /**
     * @brief Ported from corry's AnalysisEfficiency, matching the plot set
     * it actually produces (createInPixelRoiPlots() and
     * createTrackTimePlots() are commented out in corry's own source, so
     * those two groups are deliberately NOT reproduced here - only the core
     * efficiency/charge/cluster/ToA/distance plots plus the fake-rate set
     * that createFakeRatePlots() unconditionally creates).
     *
     * corry's config binds this module to detectors either via "name"
     * (exactly one detector) or "type" (auto-instantiated once per matching
     * detector - see ModuleManager.cpp) - since Warlock has one C++
     * instance per config block, "type" is handled the same way
     * DUTAssociation does: one instance internally loops over every
     * matching detector, with its own independent set of accumulators and
     * output histograms per detector (see #accum_).
     */
    class AnalysisEfficiency : public Module {
    public:
        AnalysisEfficiency(Configuration cfg, Configuration g_cfg, ThreadPool* pool);

        void initialize() override;
        void run(DataBatch& batch) override;
        void finalize() override;

    private:
        std::string dut_name_;    ///< Explicit single-detector binding ("name").
        std::string type_filter_; ///< Multi-detector binding ("type").
        std::vector<std::string> target_detectors_; ///< Resolved detector list this instance runs over - either just #dut_name_ or every detector matching #type_filter_.
        /// Subset of #target_detectors_ that split into dies
        /// (DetectorGeo::dieLabels()), mapped to their own label list -
        /// for each of these, #accum_ additionally carries a
        /// "<det_name>/<die>" entry per label alongside the combined
        /// "<det_name>" one, populated by re-running runTrackLoop()/
        /// runFakeRateLoop() with a die filter (see their `die_filter`
        /// param). A single-die board like CAEN_UZH_2 (n_pixels_y == 1) is
        /// left alone - only the combined entry exists for it.
        std::map<std::string, std::vector<std::string>> die_labels_by_det_;

        double chi2_cut_; ///< Track chi2/ndof cut applied before a track counts toward efficiency.
        std::vector<std::string> required_detectors_; ///< If non-empty, a track must have an associated cluster on every one of these to count.

        // Cut-flow parameters, matching corry's config keys/defaults.
        int masked_pixel_distance_cut_;   ///< Pixel-distance cut for rejecting tracks landing near a masked pixel.
        double spatial_cut_sensoredge_;   ///< Margin (in pixels) tracks must clear the sensor edge by to be counted.
        std::string fake_rate_method_;    ///< "radius" (default) or "edge" - which fake-hit definition runFakeRateLoop() uses.
        double fake_rate_distance_;       ///< Distance threshold used by #fake_rate_method_.
        int n_charge_bins_;               ///< Bin count for the fake-hit/cluster charge histograms.
        double charge_histo_range_;       ///< Upper edge for the fake-hit/cluster charge histograms.
        double sigma_moyal_;              ///< Sigma parameter of the Moyal-charge MPV recomputation in finalize().
        double bin_size_x_, bin_size_y_;  ///< efficiency_2D_histograms_bin_size - physical (mm) bin size for the profile maps.
        /// corry's time_cut_frameedge cut compares a track's timestamp to the
        /// start/end of its enclosing readout frame; Warlock has no
        /// per-event frame-boundary concept yet, so it's read and stored
        /// for forward-compatibility but not evaluated (at its default of
        /// 0 the cut is mathematically a no-op: `dt < 0` is never true).
        double time_cut_frameedge_;
        // XY translation / Z-rotation alignment lives in its own dedicated
        // module, AlignmentCAEN, mirroring the separation between align
        // modules (PrealignmentCAEN/AlignmentTrackChi2/AlignmentMillepede)
        // and measure modules (AnalysisEfficiency/AnalysisTiming). See
        // AlignmentCAEN's own header for the alignment logic.

        /// Per-detector accumulator bundle - one of these per entry in
        /// #target_detectors_, populated during run() (thread-locally, merged
        /// under #acc_mtx_) and turned into ratios/means in finalize().
        struct DetAccum {
            // Physical (mm-binned) efficiency maps: same bin COUNT for
            // global and local (corry quirk - both derive their bin count
            // from the same "detector_size_global" formula), global
            // centered on the detector position, local on 0.
            int n_bins_x{0}, n_bins_y{0};
            double size_x{0}, size_y{0}; // physical detector size (mm)
            Eigen::MatrixXd eff_num_global, eff_den_global;
            Eigen::MatrixXd eff_num_local, eff_den_local;

            // Local-binned "profile" accumulators (TProfile2D equivalents):
            // a shared per-bin fill count plus one sum matrix per quantity -
            // mean = sum/count at finalize(). Filled once per (track,
            // associated cluster) pair, at the track's local intercept bin
            // (not the cluster's own position) - matches corry exactly.
            Eigen::MatrixXd profile_count;
            Eigen::MatrixXd charge_sum, moyal_sum, clustersize_sum, ncols_sum, nrows_sum, toa_sum;

            // Fake-rate accumulators (createFakeRatePlots(), always active).
            Eigen::MatrixXd fake_pixel_map;
            size_t total_events{0};
            // Time-binned "profile" (fakePixelPerEventVsTime/Long): sum+count
            // per time bin, mean computed at finalize(). corry's ranges are
            // fixed (0-3000s / 0-30000s, 3000 bins each).
            std::vector<double> fake_vs_time_sum, fake_vs_time_count;
            std::vector<double> fake_vs_time_long_sum, fake_vs_time_long_count;

            // Track selection flow counters (mirrors corry's n_track/n_dut/...).
            size_t n_track{0}, n_dut{0}, n_chi2{0}, n_masked{0}, n_requirecluster{0};
            size_t total_tracks{0}, matched_tracks{0};

            // Efficiency-vs-event graph source data.
            std::map<uint64_t, std::pair<size_t, size_t>> per_event_counts;
        };
        std::map<std::string, DetAccum> accum_; ///< One DetAccum per detector in #target_detectors_.
        std::mutex acc_mtx_; ///< Guards merging each thread chunk's local accumulators into #accum_.

        /// Applies the full track selection cut-flow (DUT bounds, chi2,
        /// isolation, masked-pixel proximity, required associated clusters)
        /// and fills the efficiency/profile maps and distance-to-hit plots
        /// for every accepted track, chunked across the thread pool.
        /// `det_name` is always the underlying geometry detector (for
        /// intercept/bounds lookups); `accum_key` selects which #accum_
        /// entry (and plot dir) results are written into - either
        /// `det_name` itself (combined) or `det_name + "/" + die` for one
        /// of DetectorGeo::dieLabels(). `die_filter` is "" (every
        /// associated cluster on `det_name` counts, matching pre-split
        /// behaviour exactly) or a specific die label (only an associated
        /// cluster whose own (col, row) falls in that die, via
        /// DetectorGeo::dieOf(), counts as a match) - the track selection
        /// cuts themselves never depend on `die_filter`, only which
        /// associated clusters count as "matched" does.
        void runTrackLoop(DataBatch& batch, const std::string& det_name, const std::string& accum_key,
                           const std::string& die_filter);
        /// Fills the fake-hit-rate plots (createFakeRatePlots() in corry):
        /// per #fake_rate_method_, either clusters/pixels outside every
        /// track's widened active-area intercept ("edge") or clusters with
        /// no track intercept within #fake_rate_distance_ pixels ("radius")
        /// count as fake, chunked per event across the thread pool.
        /// `accum_key`/`die_filter` follow the same convention as
        /// runTrackLoop() - when `die_filter` is set, only clusters/pixels
        /// whose own row falls in that die are considered at all.
        void runFakeRateLoop(DataBatch& batch, const std::string& det_name, const std::string& accum_key,
                              const std::string& die_filter,
                              const std::map<uint64_t, std::vector<std::shared_ptr<Track>>>& tracks_by_event);
    };
}
#endif
