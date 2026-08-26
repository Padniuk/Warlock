# Tracking4D
**Module Type**: *GLOBAL*
**Status**: Functional

### Description
Straight-line track finding and fitting across pixel detector planes, mirroring corryvreckan's `Tracking4D`. Seeds candidate tracks from the outermost (in Z) plane pair carrying clusters in an event, extrapolates through intermediate planes to pick up matching clusters within a spatial/time cut, and refits after every match.

Can optionally dump fitted tracks and the clusters they're built from to HDF5 (`storage_file`), so downstream stages can reload already-fitted tracks via the `Reader` module instead of re-running clustering and tracking from raw hits every time. Clusters on real pixel detectors excluded from the fit (via `exclude_detectors`) are saved too, as fresh unassociated clusters ready for `DUTAssociation` to associate downstream - this includes CAEN boards, whose clusters arrive here via `ClusteringSpatial` just like any other excluded detector.

### Parameters
* `spatial_cut_abs`: `[x, y]` half-axes (mm) of the elliptical spatial-match cut. Also accepts `spatial_cut_abs_x`/`spatial_cut_abs_y` separately, each defaulting to `0.1mm`.
* `time_cut_abs`: Time cut, in ns, applied to both seed pairs and mid-plane matches. Defaults to `690000ns`.
* `min_hits_on_track`: Minimum number of clusters a candidate must carry to be accepted as a track. Defaults to `5`.
* `exclude_detectors`: Detector names and/or types left out of tracking. Accepts a single value or a list. No default (empty).
* `unique_cluster_usage`: If `true`, tracks sharing a cluster within an event are deduplicated, keeping only the best-fit one. Defaults to `false`.
* `isolated_planes` / `isolation_cut`: Parallel lists - per-detector minimum pairwise track distance (in the plane's own units) below which both tracks are marked not-isolated at that plane. Empty (default) disables isolation evaluation entirely.
* `storage_file`: If set, dumps fitted tracks (and excluded-detector clusters) to this HDF5 file for reuse via `Reader`. No default (disabled).

### Plots produced
Once per run:

* `track_chi2` / `track_chi2ndof`
* `clustersPerTrack`, `tracksPerEvent`
* `trackAngleX` / `trackAngleY`

For each detector:

* Track residual in X/Y, broken down by cluster width (1-3 pixels)
* Distance between tracks

### Usage
```toml
[Tracking4D]
min_hits_on_track = 5
spatial_cut_abs = [0.03, 0.03]
exclude_detectors = "caendt5742"
storage_file = "output/60_tracks.h5"
```
