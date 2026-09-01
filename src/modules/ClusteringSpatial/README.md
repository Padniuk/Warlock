# ClusteringSpatial
**Module Type**: *DETECTOR*
**Detector Type**: *all*
**Status**: Functional

### Description
Groups adjacent pixel hits into clusters, per detector and per event, mirroring corryvreckan's `ClusteringSpatial`. Flood-fills 8-connected neighboring pixels within each event into `Cluster` objects, computing a charge-weighted or plain geometric centroid (per `charge_weighting`) and the cluster's local/global position from the result.

For a detector split into multiple dies (`GeometryManager`'s `dieOf()`, e.g. a two-die CAEN_UZH board), the flood-fill never crosses from one die's pixels into another's, even where their `(column, row)` indices are numerically adjacent - two dies are separate physical silicon chips with a real gap between them (the geometry file's `pitch_y_gap`), so a cluster spanning both would be a grid-index artifact, not a real hit. No-op for a detector with no die split at all.

Since a CAEN waveform-digitizer channel is fed into this same clustering path as an ordinary `Pixel` (see `WaveformSelector`), adjacent-channel CAEN hits get merged into real multi-channel clusters here too, not just MIMOSA26/RD53B pixel hits.

### Parameters
* `charge_weighting`: If `true`, the cluster centroid is charge-weighted; if `false`, a plain arithmetic mean. Defaults to `true`.
* `exclude_detectors`: Detector names and/or types to leave out of clustering entirely. Accepts a single value or a list. No default (empty).
* `discard_events_missing_on`: If non-empty, an event is only clustered if every listed detector has at least one pixel in it (the intersection of their event sets). Otherwise every detector uses its own event set independently.

### Plots produced
For each detector the following plots are produced:

* `clusterSize`, `clusterSeedCharge`, `clusterCharge`
* `clusterWidthRow` / `clusterWidthColumn`
* `clusterPositionLocal` / `clusterPositionGlobal`
* `clusterTimes`, `clusterMultiplicity`
* `clusterUncertaintyX` / `clusterUncertaintyY`

### Usage
```toml
[ClusteringSpatial]
charge_weighting = true
exclude_detectors = "caendt5742"
```
