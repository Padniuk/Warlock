#ifndef WARLOCK_WAVEFORM_HPP
#define WARLOCK_WAVEFORM_HPP
#include "Pixel.hpp"
#include <vector>
#include <utility>
#include <cmath>

namespace framework {
    class Waveform : public Pixel {
    public:
        Waveform(std::string det_id, int col, int row, double ts, std::vector<double> raw_data, uint64_t ev_id = 0, int channel = -1)
            : Pixel(std::move(det_id), col, row, 1.0, ts, ev_id, channel), samples_(std::move(raw_data)) {}

        const std::vector<double>& samples() const { return samples_; }

        double t0{0.0};
        double dt{0.0};

        /// CAEN per-channel DC-offset baseline (Volts), set by
        /// EventLoaderHDF5 from the "/hits" group's
        /// "caen_dc_offset_volts_ch<N>" attribute (written by
        /// minionRawConverter - see RawConverter.cpp's caenDcOffsetVolts).
        /// #samples() holds the DRS4-corrected ADC code, not calibrated
        /// Volts (see AppendableDataset's ds_wf_samples comment); this is
        /// the one piece WaveformProcessingCAEN::run() needs to invert that
        /// same calibration back to Volts, since it can't re-parse raw-file
        /// tags itself. 0.0 (mid-scale/no-offset default) for a sensor with
        /// no such attribute, i.e. every non-CAEN waveform source.
        double dc_offset_volts{0.0};

        double baseline{0}, noise{0}, amplitude{0}, snr{0};
        double rise_time{0}, start_time{0}, integral{0};
        /// Time-over-threshold at 50% of amplitude (falling-edge 50%
        /// crossing minus rising-edge 50% crossing, both relative to
        /// baseline) and time-over-noise (same, but at the noise threshold
        /// `baseline + noise_factor*noise` instead of 50% amplitude) - see
        /// WaveformProcessingCAEN::run() for how both are found.
        double time_over_threshold{0}, time_over_noise{0};
        bool is_processed{false};

        /// (cfd fraction, start_time) pairs for every CFD value this
        /// waveform was originally processed at - populated by
        /// WaveformProcessingCAEN::run() for a freshly-computed waveform, or by
        /// the Reader module when reloading from a saved storage_file.
        /// #start_time above always mirrors whichever fraction is first in
        /// the module's configured `cfd` list (kept for backward-compatible
        /// single-value consumers); a downstream module that cares WHICH
        /// fraction it uses (e.g. WaveformSelector's own `cfd` config key)
        /// should look up its chosen value here instead via
        /// #startTimeAtCfd().
        std::vector<std::pair<double, double>> start_time_by_cfd;

        /// Looks up the start_time recorded for `cfd` (matched within a
        /// small epsilon to tolerate float round-trips through HDF5/TOML
        /// parsing), falling back to #start_time if #start_time_by_cfd was
        /// never populated (e.g. a freshly-computed waveform that hasn't
        /// gone through a save/reload cycle) or `cfd` isn't found.
        double startTimeAtCfd(double cfd) const {
            for (auto const& [c, t] : start_time_by_cfd) {
                if (std::fabs(c - cfd) < 1e-6) return t;
            }
            return start_time;
        }
    private:
        std::vector<double> samples_;
    };
}
#endif