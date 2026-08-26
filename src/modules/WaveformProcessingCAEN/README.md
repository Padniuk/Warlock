# WaveformProcessingCAEN
**Module Type**: *DETECTOR*
**Detector Type**: *CAENDT5742*
**Status**: Functional

### Description
Computes baseline, noise, amplitude, SNR, rise time, CFD start time, and charge for every waveform in a batch, mirroring corryvreckan's own `WfSignal` analysis (as used by corry's own `WaveformProcessing`). A waveform whose computed values fall outside the configured `min_`/`max_` ranges is dropped entirely - not filled into any histogram, not saved, and skipped by downstream consumers.

Can optionally break results down per (detector, channel) for configured CAEN DUT/TREF boards, and dump every waveform's computed parameters to HDF5 (`storage_file`) for downstream reuse via the `Reader` module. `cfd` accepts a list of fractions - `start_time` is recomputed and histogrammed separately for each, and all of them are saved (not just one), so a downstream `WaveformSelector` can pick which fraction it actually wants at read time.

### Parameters
* `signal_polarity`: `+1` for positive-going waveforms, `-1` for negative. Defaults to `-1`.
* `cfd`: One or more CFD fractions (e.g. `0.7` for 70%). Accepts a single value or a list.
* `sampling_interval`: Time between consecutive samples, in ns.
* `noise_factor`: Multiple of the baseline noise standard deviation used as the threshold bounding the charge integral. Defaults to `1.0`.
* `integral_time`: Charge integration window width, in ns, from the signal's rising edge. Defaults to `2.8`.
* `min_snr` / `max_snr`: SNR range to consider a valid hit. Default to `0.0` / `1e30`.
* `min_start_time` / `max_start_time`: CFD start-time range to consider a valid hit. Default to `0.0` / `1e30`.
* `min_rise_time` / `max_rise_time`: Rise-time range to consider a valid hit. Default to `0.0` / `1e30`.
* `min_integral` / `max_integral`: Charge range to consider a valid hit. Default to `0.0` / `1e30`.
* `storage_file`: If set, dumps every valid waveform's computed parameters to this HDF5 file for reuse via `Reader`. No default (disabled).

Detectors to break down per-channel are resolved automatically to every `caendt5742` detector in the geometry file.

### Plots produced
For the overall total, and again per detector and per channel, under each configured CFD fraction's own folder:

* `snr`, `amplitude`, `noise`, `baseline`, `rise_time`, `start_time`, `charge`

Flat (not per-CFD) diagnostics: `raw_snr`, `raw_amplitude`, `raw_noise`, `raw_baseline`, `raw_rise_time`, `raw_signal_start_time`, `raw_charge`.

### Usage
```toml
[WaveformProcessingCAEN]
signal_polarity = -1
cfd = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
sampling_interval = 0.2
storage_file = "output/60_waveforms.h5"
```
