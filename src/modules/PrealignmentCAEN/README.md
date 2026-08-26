# PrealignmentCAEN
**Module Type**: *DUT*
**Detector Type**: *CAENDT5742*
**Status**: Functional

### Description
X/Y prealignment for CAEN digitizer DUTs - Warlock-only, no corryvreckan counterpart. Aligns coarse-pixel CAEN boards by profiling track intercepts at each DUT's Z position and locating the half-max edges of the resulting plateau (the DUT's illuminated physical footprint), rather than correlating against a reference plane's own hit distribution.

This works because a CAEN DUT's footprint is a hard-edged rectangle set by its physical extent, not a track-correlated peak - the standard `Prealignment` module (which looks for a correlation peak) has nothing to lock onto here, but the plateau's edges can be found directly from the intercept histogram itself.

The expected physical footprint (bounding the edge search, and sanity-checking the derived width) is always computed directly from the detector's own geometry - pixel count, pitch, and die gap for a two-die board - never from a hand-typed config value that could silently drift out of sync with the `.geo` file.

### Parameters
* `dut_names`: CAEN DUTs to align, in the order they should be processed.

### Usage
```toml
[PrealignmentCAEN]
dut_names = ["CAEN_UZH_0", "CAEN_UZH_4"]
```
