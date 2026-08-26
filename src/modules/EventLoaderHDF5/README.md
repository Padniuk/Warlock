# EventLoaderHDF5
**Module Type**: *GLOBAL*
**Detector Type**: *all*
**Status**: Functional

### Description
Loads pixel hits and (optionally) waveform samples from an HDF5 file produced by `minionRawConverter`, the event source at the start of every Warlock pipeline. No corryvreckan equivalent exists - corry reads EUDAQ2 raw data directly via `EventLoaderEUDAQ2`, while Warlock reads its own pre-converted HDF5 format instead.

Builds a route table mapping each raw sensor identity to the matching detector(s) in the current geometry - a raw sensor can fan out to several Warlock detectors via a detector's `source_sensor`/`source_channels` geometry fields, e.g. one physical multi-die CAEN board split into independent detectors. Emits pixel hits always, and waveform samples too if `load_waveforms` is set and the file actually has per-hit sample data.

### Parameters
* `file_name`: Input HDF5 file(s) to load, in order. Accepts a single value or a list. Defaults to `output.h5`.
* `load_waveforms`: If `true`, also emits per-hit waveform samples (for CAEN digitizer boards) alongside pixel hits. Defaults to `false`.
* `type`: If set, only detectors of this type are loaded. Empty means every detector in the `.geo` file.
* `exclude_detectors`: Detector names and/or types left out of loading, independent of `type` - lets a stage load several types together minus one (e.g. `mimosa26`+`cmsit`, not `caendt5742`). Accepts a single value or a list.
* `number_of_events` *(global)*: Framework-wide event limit - set under `[Warlock]`, not this module's own block.
* `batch_size` *(global)*: Events per `DataBatch` - set under `[Warlock]`.

### Usage
```toml
[EventLoaderHDF5]
file_name = ["run000274.h5"]
load_waveforms = true
exclude_detectors = "caendt5742"
```
