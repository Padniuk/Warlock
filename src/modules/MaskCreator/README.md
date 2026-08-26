# MaskCreator
**Module Type**: *DETECTOR*
**Detector Type**: *all*
**Status**: Functional

### Description
Builds a per-detector pixel occupancy map across the run, then flags outlier pixels as masked, mirroring corryvreckan's `MaskCreator`. Two methods are available: `"frequency"` (a pixel is masked if its hit count exceeds `frequency_cut` times the detector's mean occupancy) and `"localdensity"` (a kernel-density comparison against each pixel's own neighborhood).

Newly-found masked pixels are appended to the mask file already referenced in the geometry file for each detector; if none is set there, a new `mask_<detector_name>.txt` is created in the output directory. Already-existing masked pixels are preserved, not overwritten. No masks are actually applied here - that happens in the event loader when reading input data.

### Parameters
* `type`: If set, only detectors of this type are masked. Empty means every pixel detector.
* `method`: `"frequency"` or `"localdensity"`. Defaults to `"frequency"`.
* `frequency_cut`: Occupancy threshold, as a multiple of the detector's mean hit rate, above which a pixel is masked. Defaults to `50.0`. Only used in `"frequency"` mode.
* `density_bandwidth`: Neighborhood half-size (in pixels) for the local density estimator. Defaults to `10`. Only used in `"localdensity"` mode.
* `sigma_above_avg_max`: Number of standard deviations above the local average at which a pixel is flagged noisy. Defaults to `5.0`. Only used in `"localdensity"` mode.
* `rate_max`: Maximum allowed per-event rate. Defaults to `1.0`. Only used in `"localdensity"` mode.
* `mask_dead_pixels`: If `true`, also masks pixels with zero recorded hits. Defaults to `false`.

### Plots produced
For each detector the following plots are produced:

* `maskmap`: map of masked pixels
* `occupancy`: 2D occupancy histogram
* `density` / `local_significance`: diagnostic plots from the `"localdensity"` computation

### Usage
```toml
[MaskCreator]
method = "frequency"
frequency_cut = 10
```
