# WaveformSelector
**Module Type**: *DETECTOR*
**Detector Type**: *CAENDT5742*
**Status**: Functional

### Description
Selects CAEN digitizer waveforms passing a set of quality cuts and turns each survivor into a `Pixel` for downstream clustering - Warlock-only, no corryvreckan counterpart (corry applies its equivalent cuts directly inside `WaveformProcessing`). Applies inclusive min/max range cuts (SNR, amplitude, noise, rise time, start time at a configured CFD fraction, charge) to every waveform in the batch.

Each accepted waveform's channel position becomes a `Pixel` fed into `ClusteringSpatial` - the same input MIMOSA26/RD53B raw hits go through - so adjacent-channel CAEN hits get merged into real multi-channel clusters rather than becoming one synthetic single-pixel cluster per waveform. One instance per `[[WaveformSelector]]` config block, optionally bound to a single detector via `name`; an unnamed instance applies its cuts to every detector present in the batch.

### Parameters
* `name`: Detector this instance applies its cuts to. Empty means every detector present in the batch.
* `cfd`: Which CFD fraction's `start_time` to use, out of whatever `WaveformProcessingCAEN` was configured with. Defaults to `0.5`.
* `min_snr` / `max_snr`: Default to `0.0` / `1e30`.
* `min_amplitude` / `max_amplitude`: Default to `0.0` / `1e30`.
* `min_noise` / `max_noise`: Default to `0.0` / `1e30`.
* `min_rise_time` / `max_rise_time`: Default to `0.0` / `1e30`.
* `min_start_time` / `max_start_time`: Applied to `start_time` at the chosen `cfd` fraction. Default to `0.0` / `1e30`.
* `min_charge` / `max_charge`: Default to `0.0` / `1e30`.
* `dut_names`: Detectors to break down per-channel in the output histograms (cuts themselves always apply to every waveform regardless of this list). Empty resolves to every `caendt5742` detector in the geometry.

### Plots produced
For the overall total, and again per configured DUT and per channel:

* `snr`, `amplitude`, `noise`, `baseline`, `rise_time`, `start_time`, `charge`

### Usage
```toml
[[WaveformSelector]]
name = "CAEN_UZH_0"
cfd = 0.5
min_snr = 18.0
min_start_time = 30.0
max_start_time = 64.0
max_rise_time = 2.0
```
