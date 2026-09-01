# DUTAssociation
**Module Type**: *DUT*
**Detector Type**: *all*
**Status**: Functional

### Description
Matches fitted tracks to clusters on one or more DUTs, mirroring corryvreckan's `DUTAssociation` module. For every track and every configured target detector, associates *every* cluster within an elliptical spatial cut and an absolute time cut - not just the closest one, matching corry's own behavior of tracking every candidate association rather than a single best match.

For a multi-die detector, candidate clusters are restricted to the die the track itself predicts (its own local intercept, not any cluster's) - Warlock-only, 2026-09-01. Without this, a track near the physical gap between two dies would be tested against both dies' boundary-row clusters, each an independent chance to pass the ellipse cut - visible as efficiency/charge maps "bleeding" further into the gap than past the true outer edge, where only one die's clusters are ever nearby. No-op for a detector with no die split.

Bound to either a single detector (`name`) or every detector of a given `type`, with one internal loop iteration per matching detector when `type` is used.

For a two-die CAEN board (a TILGAD or TREF layout in the detector geometry), every plot below is filled both combined (both dies together) and split into `/TOP` and `/BOTTOM` subpaths, so per-die association quality can be inspected independently.

### Parameters
* `name`: Single detector this instance associates against. Mutually exclusive with `type`.
* `type`: Detector type this instance associates against - one internal association loop per matching detector. Mutually exclusive with `name`.
* `spatial_cut_abs`: `[x, y]` half-axes (in mm) of the elliptical spatial cut between a track's local intercept and a candidate cluster. A single scalar is also accepted and applied to both axes.
* `time_cut_abs`: Absolute time cut (in ns) between a track and a candidate cluster's timestamp. Defaults to effectively unrestricted (`1e30`).
* `use_cluster_centre`: If `true`, compares the track intercept against the cluster centre for the spatial cut; if `false`, against the nearest pixel within the cluster. Defaults to `false`.

### Plots produced
For the DUT (and, where the board is a two-die layout, again separately for `/TOP` and `/BOTTOM`), the following plots are produced:

* Histograms of the distance in X/Y from the cluster to the pixel closest to the track, broken down by cluster width (1/2/3 pixels)
* `hCutHisto`: histogram of clusters discarded by the spatial vs. time cut
* `no_assoc_cls`: number of associated clusters per track
* 2D histograms of the local track-to-cluster distance, before cuts (`hDist_trackCluster_2D`) and after (`hDist_trackCluster_2D_assoc`), and their global-coordinate counterparts
* `hTrackPosAssoc2D` / `hTrackPosNoAssoc2D`: global track position for tracks with/without a successful association
* `hClusterPosAssoc2D`: local position of every associated cluster

### Usage
```toml
[[DUTAssociation]]
name = "CAEN_UZH_0"
spatial_cut_abs = [0.325, 0.217]
```
