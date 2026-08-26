#ifndef WARLOCK_WAVEFORMSELECTOR_HPP
#define WARLOCK_WAVEFORMSELECTOR_HPP

#include "core/Module.hpp"
#include <atomic>
#include <string>
#include <vector>

namespace framework {
    /**
     * @brief Selects CAEN digitizer waveforms passing a set of quality cuts
     * and turns each survivor into a Pixel for downstream clustering.
     *
     * Applies inclusive min_/max_ range cuts (SNR, amplitude, noise, rise
     * time, start time at a configured CFD fraction, integral,
     * time-over-threshold at 50%, time-over-noise) to every waveform in
     * DataBatch::waveforms, matching corry's
     * WaveformProcessing::is_valid_hit_() convention. Each accepted
     * waveform's channel position becomes a Pixel pushed into
     * DataBatch::pixels, feeding the exact same ClusteringSpatial path
     * MIMOSA26/RD53B raw hits go through - so adjacent-channel CAEN hits
     * get merged into real multi-channel clusters (matching corry, which
     * runs ClusteringSpatial on caendt5742 too) rather than becoming one
     * synthetic single-pixel cluster per waveform. One instance per
     * `[[WaveformSelector]]` config block, optionally bound to a single
     * detector via `name` (unnamed applies its cuts to every detector
     * present in the batch).
     */
    class WaveformSelector : public Module {
    public:
        WaveformSelector(Configuration cfg, Configuration g_cfg, ThreadPool* pool);
        void initialize() override;
        void run(DataBatch& batch) override;
        void finalize() override;

    private:
        // Target detector for this instance ("name" in the config, doubling
        // as the TOML array-of-tables instance identifier - see
        // FrameworkManager's [[WaveformSelector]] handling). Empty means
        // "apply these cuts to every detector present in batch.waveforms",
        // which is what a single un-named [WaveformSelector] block gets.
        // Multiple named instances (one per detector, each with its own
        // cuts) are configured as repeated [[WaveformSelector]] blocks.
        std::string detector_name_;

        // Inclusive range cuts, one pair per computed parameter - same
        // min_/max_ convention as corry's WaveformProcessing::is_valid_hit_().
        // Defaults are permissive (0..1e30) so an unconfigured cut never
        // rejects anything.
        double min_snr_, max_snr_;
        double min_amplitude_, max_amplitude_;
        double min_noise_, max_noise_;
        double min_rise_time_, max_rise_time_;
        double min_start_time_, max_start_time_;
        double min_integral_, max_integral_;
        double min_time_over_threshold_, max_time_over_threshold_;
        double min_time_over_noise_, max_time_over_noise_;

        // Which CFD fraction's start_time to use for the min_/max_start_time
        // cut above and for the resulting pixel's timestamp - a waveform
        // carries one start_time per fraction WaveformProcessingCAEN was
        // configured with (see Waveform::start_time_by_cfd), so this module
        // has to be told which one it actually wants (matching corry's own
        // single `cfd` config value). No implicit "first in the list"
        // fallback here on purpose - see Waveform::startTimeAtCfd().
        double cfd_;

        // Detectors to break down per-channel in the output histograms -
        // same convention as WaveformProcessingCAEN. Empty means "every
        // caendt5742 detector in the geometry" (resolved in initialize()).
        // Cuts themselves apply to every waveform regardless of this list -
        // it only controls the per-detector/per-channel histogram folders.
        std::vector<std::string> dut_names_;

        std::atomic<size_t> total_in_{0}, total_out_{0};

        bool passesCuts(const Waveform& wf) const;
        void registerChannelPlots(const std::string& dir);
        void fillChannelPlots(const std::string& dir, double snr, double amplitude, double noise,
                               double baseline, double rise_time, double start_time, double integral,
                               double time_over_threshold, double time_over_noise);
        // Turns a passing waveform into a Pixel at its channel position and
        // pushes it into batch.pixels[det_name] - the same input
        // ClusteringSpatial already consumes for MIMOSA26/RD53B raw hits,
        // so its existing neighbor-finding merges adjacent-channel CAEN
        // hits into real multi-channel clusters (matching corry, which runs
        // ClusteringSpatial on caendt5742 too) instead of one synthetic
        // single-pixel cluster per accepted waveform.
        void makePixel(DataBatch& batch, const std::string& det_name, const Waveform& wf) const;
    };
}
#endif
