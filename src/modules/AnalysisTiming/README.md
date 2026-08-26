# AnalysisTiming
**Module Type**: *GLOBAL*
**Detector Type**: *CAENDT5742*
**Status**: Functional

### Description
Performs a cross-detector timing analysis between a set of DUTs, ported from corryvreckan's own `AnalysisTiming`. For every track carrying an associated cluster on *all* configured `detectors` (and any `require_associated_cluster_on` reference detector), fills each detector's own charge-vs-ToA map plus a pairwise delta-ToA histogram for every detector pair, using each cluster's earliest-pixel timestamp.

Warlock adds a trigger-group check corry has no concept of at all: on a CAEN DT5742 board, channels 0-7 and 8-15 are sampled by two physically independent trigger groups, so comparing timing across two channels from different groups isn't physically meaningful. `require_same_trigger_group` (on by default) skips any cluster pair whose raw channels fall in different groups. Warlock also optionally supports a time-walk correction (`time_walk_method`) that corry's own module never applies (its source has a literal `// MISSING --- Time walk correction` comment).

### Parameters
* `detectors`: Names of the DUTs to correlate against each other. A track needs an associated cluster on every one to count.
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
```toml
[AnalysisTiming]
detectors = ["CAEN_UZH_4", "CAEN_UZH_2", "CAEN_UZH_3"]
require_associated_cluster_on = ["RD53B_115"]
require_same_trigger_group = true
```
