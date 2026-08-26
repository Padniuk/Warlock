# Prealignment
**Module Type**: *DETECTOR*
**Detector Type**: *all*
**Status**: Functional

### Description
Coarse translational alignment: correlates every detector's clusters against the reference plane's, then shifts each detector's position by the correlation peak/mean, per `method`. A single-pass, non-iterative first-order alignment, mirroring corryvreckan's `Prealignment` module and typically run before `AlignmentTrackChi2` or `AlignmentMillepede`. Rotational alignment is never touched.

The reference detector itself is instantiated and histogrammed like any other (matched against its own clusters) - only the actual position shift is skipped for it.

### Parameters
* `type`: If set, only detectors of this type are aligned. Empty means every pixel detector.
* `method`: `"mean"`, `"maximum"`, or `"gauss_fit"` - which statistic of the correlation histogram the shift is derived from. Defaults to `"mean"`.
* `damping_factor`: Scales the computed shift before applying it (`1.0` = full shift). Defaults to `1.0`.
* `max_correlation_rms`: RMS threshold of the correlation histogram above which a warning is logged (the shift is still applied). Defaults to `6mm`.
* `time_cut_abs`: Absolute time cut, in ns, between a cluster on the current detector and one on the reference plane. Defaults to `1e9ns`.
* `range_abs`: Correlation histogram half-range, in mm. Defaults to `10mm`.
* `fixed_planes`: Detector names excluded from shifting (still histogrammed like any other). No default (empty).

### Plots produced
For each detector the following plots are produced:

* 1D correlation histograms in X/Y against the reference plane
* 2D correlation histograms in X/Y, local and global coordinates

### Usage
```toml
[Prealignment]
max_correlation_rms = 6.0
damping_factor = 1.0
```
