# AnalysisEfficiency
**Module Type**: *DUT*
**Detector Type**: *all*
**Status**: Functional

### Description
Measures DUT tracking efficiency by comparing associated cluster positions against the interpolated track position at the DUT, ported from corryvreckan's `AnalysisEfficiency`. Efficiency is the fraction of tracks with an associated cluster over the total number of tracks intersecting the DUT (or its region-of-interest), computed in both local and global 2D binning.

Also produces charge/cluster-size/ToA profile maps and a fake-hit-rate estimate (clusters/pixels with no nearby track). For a two-die CAEN board, every accumulator is filled both combined and split into `/TOP`/`/BOTTOM`.

Warlock adds an optional Z-rotation fit (`update_z_orientation`) not present in corry: it derives each pixel's measured centroid from the passed-track distribution and fits the DUT's in-plane rotation against its nominal pixel grid, iterating `z_rotation_iterations` times since translation and rotation are coupled enough that a single pass doesn't fully converge.

### Parameters
* `name` / `type`: Bind to a single detector or every detector of a given type (mutually exclusive, same convention as `DUTAssociation`).
* `chi2ndof_cut`: Acceptance criterion for track chi2/ndof. Defaults to `3.0`.
* `require_associated_cluster_on`: Detector names which must also have an associated cluster for a track to count. Empty means no requirement.
* `masked_pixel_distance_cut`: Pixel-distance cut for rejecting tracks landing near a masked pixel. Defaults to `1`.
* `spatial_cut_sensoredge`: Margin (in pixel-pitch fractions) a track must clear the sensor edge by to be counted. Defaults to `1.0`.
* `fake_rate_method`: `"radius"` (clusters with no track within `fake_rate_distance` pixel pitches count as fake) or `"edge"` (all DUT activity in events with no track intercepting the active area, widened by `fake_rate_distance`, counts as fake). Defaults to `"radius"`.
* `fake_rate_distance`: Distance threshold (in pixel pitches) used by `fake_rate_method`. Defaults to `2.0`.
* `n_charge_bins` / `charge_histo_range`: Bin count and upper edge for the charge histograms. Default to `1000` / `1000.0`.
* `efficiency_2D_histograms_bin_size`: `[x, y]` physical bin size (mm) for the local/global efficiency maps. Defaults to `10um` on both axes.
* `update_z_orientation`: Warlock-only. If `true`, fits and applies each DUT's Z-rotation correction, writing the result to `detectors_file_updated`. Defaults to `false`.
* `z_rotation_iterations`: Upper bound on rebin+refit passes `update_z_orientation` performs. Defaults to `12`.
* `min_quadrant_entries`: Minimum weighted passed-histogram entries a pixel needs before its centroid is trusted by the rotation fit. Defaults to `50.0`.
* `z_rotation_row`: Restricts the rotation fit to one physical die - `0` uses every pixel, positive uses only the top row, negative only the bottom row. Defaults to `0`. When left at `0`, a dead/near-empty pixel is auto-detected (its fitted centroid lands more than half a pitch from its own nominal position) and the fit auto-restricts to the other row instead - an explicit non-zero value here always overrides that.

### Plots produced
For each detector (and, where split, again for `/TOP`/`/BOTTOM`):

* `globalEfficiencyMap_trackPos` / `localEfficiencyMap_trackPos`: efficiency as a ROOT-style passed/total ratio, global and local coordinates
* `localChargeMapMean_trackPos`, `localMoyalCharge_trackPos`, `localMPVCharge_trackPos`
* `clusterSize_trackPos`, `clusterNCols_trackPos`, `clusterNRows_trackPos`, `ToAMap_trackPos`
* `fake_rate/fakePixelPerEventMap` and its time-binned profiles
* `efficiency_vs_event`, `distanceTrackHit`, `distanceTrackHit2D`

### Usage
```toml
[AnalysisEfficiency]
type = "caendt5742"
chi2ndof_cut = 6
require_associated_cluster_on = "RD53B_115"
```
