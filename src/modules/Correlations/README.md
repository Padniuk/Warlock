# Correlations
**Module Type**: *DETECTOR*
**Detector Type**: *all*
**Status**: Functional

### Description
Fills raw hitmaps and pixel/cluster correlation histograms (position and time, local and global) between every detector and the reference plane, mirroring corryvreckan's `Correlations` module - the standard first diagnostic run after loading data, used to check gross alignment/timing before any actual alignment module runs.

The reference detector itself is included and produces real self-correlation output (against its own hits/clusters), matching corry's own behavior.

### Parameters
* `do_time_cut`: If `true`, cluster pairs outside `time_cut_abs` are excluded from the position-correlation fills. Defaults to `false`.
* `time_cut_abs`: Absolute cluster-pair time cut, in ns, used only when `do_time_cut` is set. Defaults to `690000ns`.
* `time_binning`: Bin width of the correlation-time histogram, in ns. Defaults to `1.0ns`.

### Plots produced
For each device the following plots are produced:

* 2D: hitmaps (pixel and cluster level), time correlation over time and over raw value, spatial correlations in X/Y/columns/rows (local coordinates), X/Y correlations (global coordinates)
* 1D: X/X, Y/Y, X/Y, Y/X correlations against the reference, correlation-time histograms (pixel and cluster level)

### Usage
```toml
[Correlations]
do_time_cut = true
time_cut_abs = 200ns
```
