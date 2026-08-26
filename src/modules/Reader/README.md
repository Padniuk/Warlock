# Reader
**Module Type**: *GLOBAL*
**Status**: Functional

### Description
Reloads tracks, clusters, and waveform parameters previously saved by `Tracking4D`/`WaveformProcessingCAEN` back into a `DataBatch` - Warlock-only, no corryvreckan counterpart. Lets a downstream analysis stage (waveform selection, DUT association, timing) reuse already-fitted/-clustered data instead of re-running clustering and tracking from raw hits every time.

Reads whichever of `/tracks`, `/track_cluster_links`, `/clusters/<name>`, `/waveforms/<name>` each input file actually contains. Clusters never referenced by a track link (e.g. a DUT excluded from the fit) are replayed as fresh, unassociated clusters, ready for `DUTAssociation` to associate downstream. All rows are loaded eagerly at startup; `run()` then replays them batch by batch, one track-event-range at a time.

### Parameters
* `filenames`: Input HDF5 file(s) to load. Accepts a single value or a list.
* `batch_size` *(global)*: Tracks emitted per `run()` call, capped by this - set under `[Warlock]`, not this module's own block.

### Usage
```toml
[Reader]
filenames = ["output/60_tracks.h5"]
```
