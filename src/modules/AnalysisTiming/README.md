# AnalysisTiming
**Module Type**: *GLOBAL*
**Detector Type**: *CAENDT5742*
**Status**: Functional

### Description
Performs a cross-detector timing analysis between a set of DUTs, ported from corryvreckan's own `AnalysisTiming`. Every track carrying an associated cluster on any `require_associated_cluster_on` reference detector is considered; for each detector pair in `detectors`, that pair's own charge-vs-ToA map and pairwise delta-ToA histogram are filled from tracks meeting `require_cluster_on_all_detectors`' condition (see below) - by default, a track just needs a cluster on *both* detectors of that specific pair, not on every configured detector (Warlock-only deviation from corry, 2026-08-31: with `detectors` split across disjoint DT5742 trigger groups, requiring a cluster on every single one starves every pair within the less-populated group of nearly all its statistics, since very few tracks physically cross every configured board's less-populated die at once).

Warlock adds a trigger-group check corry has no concept of at all: on a CAEN DT5742 board, channels 0-7 and 8-15 are sampled by two physically independent trigger groups, so comparing timing across two channels from different groups isn't physically meaningful. `require_same_trigger_group` (on by default) skips any cluster pair whose raw channels fall in different groups. Warlock also optionally supports a time-walk correction (`time_walk_method`) that corry's own module never applies (its source has a literal `// MISSING --- Time walk correction` comment).

### Parameters
* `detectors`: Names of the DUTs to correlate against each other. Whether a track needs a cluster on every one of these (vs. just each pair it contributes to) is controlled by `require_cluster_on_all_detectors`.
* `require_cluster_on_all_detectors`: Warlock-only. `false` (default): each detector pair is gated independently - a track only needs a cluster on both detectors of that specific pair. `true`: a track needs a cluster on *every* detector in `detectors` before it counts toward any pair at all. Needed for a proper 3-detector ("triplet") timing-resolution measurement - solving `sigma_i^2 + sigma_j^2 = sigma_ij^2` for each detector's own resolution requires all three of a triplet's pairwise dtoa distributions to come from the same track population, not three independently-selected ones. The intended pattern: one `[[AnalysisTiming]]` instance per triplet, `detectors` set to exactly those 3, `require_cluster_on_all_detectors = true`, each instance given its own `name` so their `dtoa_*` histograms land in separate output paths instead of colliding.
* `require_associated_cluster_on`: Additional detector(s) a track must also have an associated cluster on, without being correlated themselves (e.g. a reference plane used only as a gate).
* `require_same_trigger_group`: Warlock-only. If `true`, a cluster pair is only compared when both clusters' raw CAEN channels fall in the same DT5742 trigger group. Defaults to `true`.
* `time_walk_method`: `"none"` (default), `"inverse_charge"`, or `"inverse_sqrt_charge"` - Warlock-only, off by default, corry never applies this.
* `time_walk_a_<detector>` / `time_walk_b_<detector>`: Per-detector time-walk coefficients, used only when `time_walk_method` is set. Default to `0.0`.

### Plots produced
For each detector (and again per `/TOP`/`/BOTTOM` die where applicable):

* `charge_toa_<det>`: cluster charge vs. raw time-of-arrival
* `npairs_<det>`: number of cluster pairs from this detector associated to the same track
* `ncat_<det>`: number of clusters from this detector associated to the same track

For each detector pair:

* `dtoa_<det1>_vs_<det2>`: difference in time-of-arrival
* `timewalk_<det1>_in_<det1>_vs_<det2>`: amplitude vs. delta-ToA correlation (Warlock-only, useful for deriving time-walk parameters offline)

### Usage
Plain per-pair mode (default):
```toml
[AnalysisTiming]
detectors = ["CAEN_UZH_4", "CAEN_UZH_2", "CAEN_UZH_3"]
require_associated_cluster_on = ["RD53B_115"]
require_same_trigger_group = true
```

Triplet mode - one instance per triplet, each independently gated:
```toml
[[AnalysisTiming]]
name = "triplet_012"
detectors = ["CAEN_UZH_0", "CAEN_UZH_1", "CAEN_UZH_2"]
require_cluster_on_all_detectors = true
require_associated_cluster_on = ["RD53B_115"]
require_same_trigger_group = true

[[AnalysisTiming]]
name = "triplet_013"
detectors = ["CAEN_UZH_0", "CAEN_UZH_1", "CAEN_UZH_3"]
require_cluster_on_all_detectors = true
require_associated_cluster_on = ["RD53B_115"]
require_same_trigger_group = true
```
