# The Warlock Grimoire

```
   _      __         __         __
  | | /| / /__ _____/ /__  ____/ /__
  | |/ |/ / _ `/ __/ / _ \/ __/  '_/
  |__/|__/\_,_/_/ /_/\___/\__/_/\_\
  Reconstruction Framework
```

This is the reference documentation for **Warlock**, a from-scratch C++
reimplementation of a corryvreckan-style particle-track reconstruction
framework. It reads raw pixel/waveform hit data, clusters it, fits tracks
through a telescope of detector planes, and produces the same class of
physics output as the reference framework it was built to match.

## Reading order

A run is one @ref framework::FrameworkManager, which parses a TOML config,
constructs a pipeline of @ref framework::Module instances (one per
`[ModuleName]` block, in the order they're declared), and then repeatedly:

1. Loads one batch of events into a @ref framework::DataBatch (bounded by
   `batch_size`, since a full run's data doesn't fit in memory at once).
2. Runs every module's `run(DataBatch&)` over that batch, in pipeline order
   - each module reads fields another upstream module wrote and writes
   fields a downstream one will read (pixels → clusters → tracks; waveforms
   are the one independent branch, see @ref framework::Module::runsConcurrently).
3. Repeats until the configured event/track limit is reached.

Everything CPU-heavy runs across a single shared @ref framework::ThreadPool,
constructed once for the whole run and handed to every module at
construction time.

## Where things live

- @ref framework::Module / @ref framework::DataBatch - the pipeline contract
  every reconstruction module is built against.
- @ref framework::FrameworkManager - owns the pipeline, the thread pool, and
  the batch loop.
- @ref framework::Configuration - typed access to a module's own TOML block.
- @ref framework::GeometryManager / `Detector` - the telescope geometry
  (plane positions, pixel pitch, masked pixels) every module queries.
- @ref framework::PlotManager - the single histogram/graph registry every
  module fills into; written out to HDF5 at the end of a run (and,
  optionally, converted to ROOT afterward by the separate `minionRootPlotter`
  tool - Warlock itself never links against ROOT).
- @ref framework::AppendableDataset - chunked, incrementally-growable HDF5
  columns, used anywhere a module streams per-event rows to disk instead of
  holding them all in memory.
- @ref framework::Logger / `WR_LOG` - the framework's own leveled, aligned,
  colored logging.

## A concrete module

@ref framework::MaskCreator is a good first module to read end to end: it
takes raw pixel hits, builds an occupancy map per detector, flags
statistical outlier pixels as masked, and writes both back out. It touches
most of the same machinery (`Configuration`, `GeometryManager`,
`PlotManager`, HDF5 output) that every other module does, without the added
complexity of multithreaded fitting or waveform math.
