/**
 * @file Synchronizer.hpp
 * @brief Detects (and, once you've decided, corrects) a CAEN/MIMOSA
 * event-index desync directly from already-computed Tracking4D and
 * WaveformProcessingCAEN output.
 *
 * Three-stage design:
 *
 * 1. CUSUM changepoint detection, run independently per detector, finds
 *    each detector's own candidate onset. A single detector's noise
 *    characteristics can make this fire too early, so one detector alone
 *    isn't trustworthy for exactly where the break occurred.
 *
 * 2. Combine to a shared onset: every detector in a group shares one
 *    physical digitizer, so a real desync must hit all of them at the
 *    same true event. The group's onset is the max (not min or average)
 *    of the individually-detected break points, since noise can only
 *    inflate CUSUM's cumulative sum, never deflate a real break below its
 *    guaranteed detection - a too-early trigger is far more likely than a
 *    too-late one. Detectors that never cross CUSUM's threshold are
 *    excluded from this max entirely: "never fired" means insufficient
 *    evidence, not confirmed clean.
 *
 * 3. From that shared onset, a banded sequential shift search pools
 *    hit/total counts across every detector in the group for each
 *    candidate shift, rather than trusting any single detector's noisy
 *    SNR-based signal alone. The search band re-centers on each window's
 *    own answer (a Sakoe-Chiba band, as in Dynamic Time Warping), since
 *    the true shift can grow across the file rather than settling to one
 *    constant.
 *
 * The per-event pass/fail signal replicates WaveformSelector's cuts and
 * DetectorGeo::hasIntercept()'s acceptance check, so detection doesn't
 * depend on AnalysisEfficiency having been run - only
 * WaveformProcessingCAEN's storage_file and Tracking4D's save_tracks_file
 * are needed.
 */
#ifndef WARLOCK_SYNCHRONIZER_HPP
#define WARLOCK_SYNCHRONIZER_HPP

#include <cstdint>
#include <string>
#include <vector>

namespace framework {

    struct SynchronizerConfig {
        std::string tracks_file;      ///< Tracking4D's save_tracks_file (MIMOSA-only fitted tracks).
        std::string waveforms_file;   ///< WaveformProcessingCAEN's storage_file.
        std::string geo_file;         ///< The .toml geometry file the run used - real acceptance geometry.

        /// The detector GROUP to analyze jointly - every detector sharing
        /// one physical digitizer, e.g. every CAEN_IJS_* board. A single
        /// entry behaves exactly like independent single-detector
        /// analysis (the shared-onset max and the pooled shift search
        /// both reduce to the single-detector case with one member).
        std::vector<std::string> detectors;

        double spatial_cut_sensoredge = 0.0;

        // WaveformSelector-equivalent cuts - same names, same defaults
        // (no-op bounds), same semantics as WaveformSelector.cpp's own
        // config options, applied identically to every detector in the
        // group (reasonable since a shared digitizer's boards normally
        // share one selection config too - pass different Synchronizer
        // instances if that's ever not true).
        double min_snr = 0.0, max_snr = 1e30;
        double min_amplitude = 0.0, max_amplitude = 1e30;
        double min_noise = 0.0, max_noise = 1e30;
        double min_rise_time = 0.0, max_rise_time = 1e30;
        double min_start_time = 0.0, max_start_time = 1e30;
        double min_charge = 0.0, max_charge = 1e30;
        double cfd = 0.5;  ///< Must match one of storage_file's saved cfd_fractions exactly.

        // CUSUM changepoint detection (per detector).
        size_t cusum_baseline_tracks = 5000;
        double cusum_slack = 0.01;
        double cusum_threshold = 15.0;

        // Banded sequential shift search (pooled across the group, from
        // the shared onset onward).
        size_t shift_window_tracks = 5000;  ///< Pooled in-acceptance tracks per window, summed across the group.
        int64_t shift_band = 30;            ///< Search +/- this many shifts around the PREVIOUS window's answer.
    };

    /// One window's own best-fit shift, pooled across every detector in
    /// the group, covering the closed event_id range [event_lo, event_hi].
    struct ShiftSegment {
        uint64_t event_lo = 0, event_hi = 0;
        int64_t shift = 0;
        double recovered_rate = 0.0;   ///< Pooled pass rate at `shift` within this window.
        double unshifted_rate = 0.0;   ///< Pooled pass rate at shift=0 within this window, for comparison.
    };

    /// One detector's own independent CUSUM result - stage 1 only, before
    /// any cross-detector combination.
    struct PerDetectorBreak {
        std::string detector;
        bool shift_tested = false;   ///< False if geometry/waveforms_file didn't have usable data for it at all.
        bool found_break = false;    ///< True only if CUSUM actually crossed threshold - NOT evidence of "clean" when false.
        double baseline_pass_rate = 0.0;
        uint64_t break_event = 0;
    };

    /// The group's combined result: every member's own independent
    /// PerDetectorBreak, the shared onset derived from them, and the
    /// pooled shift trajectory from that onset onward.
    struct GroupResult {
        std::vector<PerDetectorBreak> per_detector;
        bool found_break = false;          ///< True if at least one member found a break.
        uint64_t shared_break_event = 0;   ///< max(break_event) over members with found_break == true.
        std::vector<ShiftSegment> segments;
    };

    class Synchronizer {
    public:
        explicit Synchronizer(SynchronizerConfig cfg);

        /// Read-only: runs all three stages for cfg_.detectors as one
        /// group. Never opens waveforms_file for writing.
        GroupResult detect();

        /// Every detector present in both the loaded geometry and
        /// waveforms_file - for CLI use when the caller wants to discover
        /// candidates rather than name a group explicitly. Requires
        /// detect()/loadGeometry to have already run, or call
        /// loadGeometryForDiscovery() first.
        std::vector<std::string> availableDetectors() const;
        void loadGeometryForDiscovery() const;

        /// Applies a piecewise shift schedule to EVERY detector in
        /// `detectors`' event_id column in waveforms_file (the same
        /// segments for all of them - that's the point of a shared
        /// group): a row with event_id inside some segment's
        /// [event_lo, event_hi] gets that segment's shift added; rows
        /// outside every segment are left untouched. Writes the WHOLE
        /// file (every detector present, not just the group) to
        /// output_file. Refuses if output_file == waveforms_file unless
        /// allow_in_place is true.
        void applyShift(const std::vector<std::string>& detectors, const std::vector<ShiftSegment>& segments,
                         const std::string& output_file, bool allow_in_place);

    private:
        SynchronizerConfig cfg_;

        struct TrackRow {
            uint64_t event_id;
            double mx, cx, my, cy;
        };
        std::vector<TrackRow> loadTracks() const;

        struct WfRow {
            uint64_t event_id;
            double snr, amplitude, noise, rise_time, start_time, charge;
        };
        std::vector<WfRow> loadWaveforms(const std::string& detector) const;
        bool passesCuts(const WfRow& w) const;
    };
}
#endif
