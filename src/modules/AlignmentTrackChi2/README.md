# AlignmentTrackChi2
**Module Type**: *GLOBAL*
**Status**: Functional

### Description
Iterative translational and rotational telescope-plane alignment, architecturally similar to corryvreckan's own `AlignmentTrackChi2`: for each non-reference, non-fixed detector in turn, an Eigen Levenberg-Marquardt fit perturbs that detector's position/orientation and refits every track touching it, minimizing `sqrt(chi2)` per track. This repeats for `iterations` full passes over all detectors.

Warlock uses a different numerical optimizer than corry's own implementation, so convergence is not bit-exact between the two - agreement should be judged by the final detector position, not by intermediate values.

### Parameters
* `iterations`: Number of full passes over all detectors. Defaults to `5`.
* `align_position`: Whether the detector's X/Y position is a free alignment parameter (Z is never aligned). Defaults to `true`.
* `align_orientation`: Whether the detector's X/Y/Z rotation is a free alignment parameter. Defaults to `true`.
* `prune_tracks`: If `true`, tracks whose chi2/ndof exceeds a fixed internal threshold are excluded from the fit. Defaults to `false`.
* `type`: If set, only detectors of this type are aligned. Empty means every detector.
* `fixed_planes`: Detector names excluded from alignment despite not being the reference plane. No default (empty).
* `number_of_tracks` *(global)*: Framework-wide track limit, used only to report progress - set under `[Warlock]`, not this module's own block.

### Plots produced
For each aligned detector, the following plots are produced:

* Graphs of the translational shift along X/Y vs. iteration number
* Graphs of the rotational shift along X/Y/Z vs. iteration number

### Usage
```toml
[AlignmentTrackChi2]
iterations = 5
align_position = true
align_orientation = true
```
