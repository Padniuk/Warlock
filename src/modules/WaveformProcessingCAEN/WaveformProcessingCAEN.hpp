/**
 * @file WaveformProcessingCAEN.hpp
 * @brief Per-waveform baseline/noise/amplitude/CFD-timing extraction for digitizer channels.
 */

#ifndef WARLOCK_WAVEFORMPROCESSINGCAEN_HPP
#define WARLOCK_WAVEFORMPROCESSINGCAEN_HPP

#include "core/Module.hpp"
#include "core/utils/AppendableDataset.hpp"
#include <H5Cpp.h>
#include <atomic>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace framework {
    class Histogram1D;
    /**
     * @brief Computes baseline, noise, amplitude, SNR, rise time, CFD start
     * time, integral, and time-over-threshold/time-over-noise for every
     * waveform in a batch, mirroring corryvreckan's WfSignal analysis for
     * everything but the latter two (Warlock-only additions - see run() for
     * the per-quantity porting notes).
     * Chunked across the thread pool per detector (sequential relative to
     * the rest of the pipeline, though - see runsConcurrently()'s docs for
     * why: WaveformSelector runs downstream in the same pipeline and
     * mutates DataBatch::waveforms too). Optionally breaks results down per
     * (detector, pixel) channel for configured DUT/TREF boards, and can dump
     * every waveform's computed parameters to HDF5 for downstream reuse
     * (see WaveformSelector / the Reader module).
     */
    class WaveformProcessingCAEN : public Module {
    public:
        WaveformProcessingCAEN(Configuration cfg, Configuration g_cfg, ThreadPool* pool);
        void initialize() override;
        void run(DataBatch& batch) override;
        void finalize() override;

        // No override here on purpose: runsConcurrently() (see Module.hpp)
        // reads the `run_concurrently` config key, defaulting to false. Do
        // NOT set it true in any config that also runs WaveformSelector
        // downstream in the same pipeline: that module both reads and
        // mutates batch.waveforms (`waveforms = std::move(accepted)`), and
        // FrameworkManager launches a runsConcurrently() module on its own
        // detached thread without waiting - so WaveformSelector could run
        // alongside this module's still-in-flight worker threads, causing
        // a use-after-free on the moved-from waveform vector. Safe to set
        // true only in a config where nothing downstream touches
        // batch.waveforms again.

    private:
        int polarity_; ///< Signal polarity (+1 or -1); samples are multiplied by this before analysis so the pulse always points "up".
        std::vector<double> cfd_values_; ///< One or more CFD fractions - start_time is recomputed and histogrammed separately for each. A single scalar `cfd = 0.7` in the config is treated as a one-element list.
        double sampling_interval_; ///< Time between consecutive samples, in ns.
        double noise_factor_; ///< Threshold (as a multiple of the baseline noise std-dev) used to bound the charge integral - see run() for the exact formula.
        double integral_time_; ///< Charge integration window width (ns) from the signal's rising edge, matching corry's `integral_time` (default 2.8, see WaveformProcessingCAEN::getSignalIntegralDt()).

        /// Hit-selection thresholds, ported from corry's `is_valid_hit_()` -
        /// a waveform whose computed snr/start_time/rise_time/integral fall
        /// outside these ranges is dropped entirely (wf->is_processed stays
        /// false: not filled into any histogram, not saved, and skipped by
        /// downstream consumers like WaveformSelector - see run()). Config
        /// keys and defaults match corry's exactly (`min_snr`, `max_snr`,
        /// `min_start_time`, `max_start_time`, `min_rise_time`,
        /// `max_rise_time`, `min_integral`, `max_integral`), so a config
        /// that doesn't set any of them reproduces corry's own default
        /// behavior: reject only if snr/start_time/rise_time/integral end
        /// up at exactly 0 (e.g. no valid signal found), not real cuts.
        double min_snr_, max_snr_, min_start_time_, max_start_time_,
               min_rise_time_, max_rise_time_, min_integral_, max_integral_;
        /// Detectors to break down per-channel (e.g. every CAEN DUT/TREF
        /// board) - each gets its own "<det_name>/total" (all its channels
        /// combined) plus one "<det_name>/channel_<n>" folder per physical
        /// pixel, nested under every configured CFD value's own folder.
        /// Detectors not listed here still contribute to the overall "total".
        /// Empty means "every caendt5742 detector in the geometry file"
        /// (resolved in initialize()).
        std::vector<std::string> dut_names_;
        std::atomic<size_t> total_processed_{0}; ///< Cumulative waveform count across every batch, reported in finalize().

        static std::string cfdDirName(double cfd);
        /// Registers the same fixed set of histograms (snr, amplitude, noise,
        /// baseline, rise_time, start_time, integral, time_over_threshold,
        /// time_over_noise) under one path - called once for the run-wide
        /// total, once per DUT's total, once per DUT per channel, for every
        /// configured CFD value.
        void registerChannelPlots(const std::string& dir);

        /// Direct pointers to one channel's 9 histograms (snr, amplitude,
        /// noise, baseline, rise_time, start_time, integral,
        /// time_over_threshold, time_over_noise), resolved once in
        /// initialize() right after registerChannelPlots(dir) - so run()'s
        /// per-waveform, per-CFD-value hot loop can fill them directly
        /// instead of rebuilding "<dir>/<name>" strings and re-doing a
        /// PlotManager map lookup per histogram.
        /// All members stay null together if `dir` was never registered
        /// (e.g. a configured dut_names_ entry absent from the geometry
        /// file) - fillChannelPlots() checks the first pointer and no-ops in
        /// that case, matching PlotManager::fill1D()'s own
        /// "unregistered -> no-op" contract.
        struct ChannelPlotRefs {
            Histogram1D *snr = nullptr, *amplitude = nullptr, *noise = nullptr, *baseline = nullptr,
                        *rise_time = nullptr, *start_time = nullptr, *integral = nullptr,
                        *time_over_threshold = nullptr, *time_over_noise = nullptr;
        };
        ChannelPlotRefs fetchChannelPlots(const std::string& dir);
        void fillChannelPlots(const ChannelPlotRefs& refs, double snr, double amplitude, double noise,
                               double baseline, double rise_time, double start_time, double integral,
                               double time_over_threshold, double time_over_noise);
        /// Sets `title` on every histogram in `refs` (see
        /// Histogram1D::setTitle()) - used to stamp each per-channel
        /// folder's plots with which raw CAEN hardware channel they
        /// actually are (channel_<col + n_pixels_x*row>, the pixel-grid
        /// index used for the folder name, is NOT the same number as the
        /// physical CH<n> label printed in RawConverter's own startup log
        /// - e.g. for CAEN_UZH_4/TREF, channel_0 is raw CH8 and channel_1
        /// is raw CH0). Only #fillChannelPlots's callers know the raw
        /// channel (it's per-waveform data, not known at registration
        /// time), so this gets called lazily on first fill rather than
        /// during initialize() - see #dut_channel_titled_.
        void titleChannelPlots(const ChannelPlotRefs& refs, const std::string& title);

        /// dut_names_[i] -> i, resolved once in initialize() so run() can
        /// find a waveform's DUT slot (into #dut_total_refs_ /
        /// #dut_channel_refs_) in O(1) instead of re-scanning dut_names_
        /// per waveform.
        std::unordered_map<std::string, size_t> dut_index_;
        std::vector<ChannelPlotRefs> cfd_total_refs_;                        ///< [cfd_idx]
        std::vector<std::vector<ChannelPlotRefs>> dut_total_refs_;           ///< [dut_idx][cfd_idx]
        std::vector<std::vector<std::vector<ChannelPlotRefs>>> dut_channel_refs_; ///< [dut_idx][cfd_idx][channel]
        /// Parallel to #dut_channel_refs_ - true once that channel's plots
        /// have been titled with their real raw hardware channel (see
        /// #titleChannelPlots()). Plain uint8_t, not atomic/mutex-guarded:
        /// worst case under a race is a few redundant (but identical, so
        /// harmless) setTitle() calls from different threads before every
        /// worker observes true - never a torn read, since a byte write is
        /// already atomic on every platform this runs on, and the exact
        /// moment it's observed doesn't affect correctness, only how many
        /// redundant calls happen right at the start of a run.
        std::vector<std::vector<std::vector<uint8_t>>> dut_channel_titled_; ///< [dut_idx][cfd_idx][channel]

        /// Optional: dump per-waveform computed parameters (baseline, noise,
        /// amplitude, snr, rise_time, start_time, integral,
        /// time_over_threshold, time_over_noise) to an HDF5 file, keyed by
        /// (event_id, col, row, channel) per detector, so
        /// downstream stages (WaveformSelector, ...) can reload already-
        /// computed parameters via the Reader module instead of recomputing
        /// them from raw samples every time. `channel` is the raw hardware
        /// channel (see Pixel::channel()), -1 if unavailable - lets a
        /// downstream module distinguish which physical die of a multi-pad
        /// board a saved waveform came from without needing separate
        /// detector geometry per die.
        std::string storage_file_;
        bool save_enabled_{false};
        std::unique_ptr<H5::H5File> save_h5_;

        /// Per-detector set of AppendableDatasets backing one
        /// "/waveforms/<det_name>" HDF5 group.
        struct ChannelWriter {
            std::unique_ptr<AppendableDataset> event_id, col, row, channel, timestamp,
                baseline, noise, amplitude, snr, rise_time, integral, t0, dt,
                time_over_threshold, time_over_noise;
            /// Rank-2 dataset, one column per configured CFD fraction (see
            /// #cfd_values_) - width recorded via the "cfd_values" HDF5
            /// attribute on the same group so a reader can match columns
            /// back to fractions without assuming an order/count.
            std::unique_ptr<AppendableDataset> start_time_cfd;
        };
        std::map<std::string, ChannelWriter> channel_writers_;

        /// Lazily creates (or returns the existing) ChannelWriter for `det_name`.
        ChannelWriter& getChannelWriter(const std::string& det_name);
    };
}
#endif
