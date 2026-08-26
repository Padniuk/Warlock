# AlignmentMillepede
**Module Type**: *GLOBAL*
**Status**: Functional

### Description
Global track-based alignment via the classic Millepede-II algorithm, ported from corryvreckan's own `AlignmentMillepede` (originally from the [Kepler framework](https://gitlab.cern.ch/lhcb/Kepler) used by LHCb). Every track contributes local (track-parameter) and global (alignment-parameter) measurement equations; local parameters are eliminated per track, and the resulting reduced global normal equations are solved once per outer iteration to get a simultaneous correction for every plane at once - unlike `AlignmentTrackChi2`, which optimizes one detector at a time.

The linear-algebra routines are direct ports of Millepede-II's own Gauss-Jordan elimination, kept structurally close to the original rather than replaced with a library solver so intermediate behavior (including corry's own quirks) matches as closely as possible.

The module stops iterating once convergence - the mean absolute correction across all free parameters - drops below `convergence`.

### Parameters
* `exclude_detectors`: Detector names and/or types to leave out of the alignment. Accepts a single value or a list. No default (empty).
* `iterations`: Number of outer (global) Millepede iterations. Defaults to `5`.
* `dofs`: Six boolean flags selecting which alignment degrees of freedom are free - Translation X/Y/Z, Rotation X/Y/Z, in that order. Defaults to `true, true, false, true, true, true` (Z translation held fixed).
* `sigmas`: Six prior uncertainties for the DOFs above, in their respective units. Defaults to `50um, 50um, 500um, 5mrad, 5mrad, 5mrad`.
* `residual_cut`: Per-equation residual cut used from the second inner iteration onward. Defaults to `0.05mm`.
* `residual_cut_init`: Looser residual cut applied only on the very first inner iteration. Defaults to `0.6mm`.
* `number_of_stddev`: If nonzero, tracks whose local chi2 exceeds this many standard deviations are rejected as outliers. Defaults to `0` (disabled).
* `convergence`: Outer-iteration loop stops once the mean absolute parameter correction drops below this value. Defaults to `1e-5`.

### Usage
```toml
[AlignmentMillepede]
iterations = 10
dofs = [true, true, false, true, true, true]
exclude_detectors = "RD53B_115"
```
