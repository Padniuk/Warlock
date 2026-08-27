# Reader
**Module Type**: *GLOBAL*
**Status**: Functional

### Description
Reloads tracks, clusters, and waveform parameters previously saved by `Tracking4D`/`WaveformProcessingCAEN` back into a `DataBatch` - Warlock-only, no corryvreckan counterpart. Lets a downstream analysis stage (waveform selection, DUT association, timing) reuse already-fitted/-clustered data instead of re-running clustering and tracking from raw hits every time.

Reads whichever of `/tracks`, `/track_cluster_links`, `/clusters/<name>`, `/waveforms/<name>` each input file actually contains. Clusters never referenced by a track link (e.g. a DUT excluded from the fit) are replayed as fresh, unassociated clusters, ready for `DUTAssociation` to associate downstream. All rows are loaded eagerly at startup; `run()` then replays them batch by batch, one track-event-range at a time.

### Parameters
* `filenames`: Input HDF5 file(s) to load. Accepts a single value or a list.
* `max_event_id` *(optional)*: Rows with `event_id` beyond this are silently dropped at load time - stop replaying this run at a known-good point (e.g. right before a detected desync) without discarding everything before it and without re-running the upstream tracking stage that produced the file. Unlimited by default. Local to this `[[Reader]]` block, not `[Warlock]`'s global `number_of_events` - each instance is one independent run with its own cutoff.
* `batch_size` *(global)*: Tracks emitted per `run()` call, capped by this - set under `[Warlock]`, not this module's own block.

### Usage
```toml
[Reader]
filenames = ["output/60_tracks.h5"]
```
