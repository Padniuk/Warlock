/**
 * @file Module.hpp
 * @brief The reconstruction-module base class and its shared per-batch state.
 */

#ifndef WARLOCK_MODULE_HPP
#define WARLOCK_MODULE_HPP

#include "core/config/Configuration.hpp"
#include "core/utils/ThreadPool.hpp"
#include "objects/Pixel.hpp"
#include "objects/Cluster.hpp"
#include "objects/Track.hpp"
#include "objects/Waveform.hpp"

#include <map>
#include <vector>
#include <memory>
#include <string>
#include <filesystem>

namespace framework {

    /**
     * @brief Per-batch working set shared by every module in the pipeline.
     *
     * One instance is constructed per batch of events and passed by
     * reference through every module's run() in pipeline order. Modules
     * communicate purely through the fields below: an upstream module
     * writes a field, a downstream one reads it (e.g. an event-source
     * module fills #pixels, a clustering module reads #pixels and fills
     * #clusters, a tracking module reads #clusters and fills #tracks).
     */
    struct DataBatch {
        /// Running count of events processed before this batch; also used
        /// as this batch's own starting event id.
        size_t batch_id;

        /// Set by a module to request the pipeline stop after this batch.
        bool halt_pipeline{false};

        /// Requested size (in events) for the *next* batch, or 0 to use the
        /// configured `batch_size` as-is. Set by FrameworkManager as a
        /// `number_of_tracks` target gets close, so the final batch doesn't
        /// read/cluster/track a full extra chunk of events beyond what's
        /// actually needed to reach it.
        size_t requested_batch_size{0};

        /// Actual number of distinct events the event-source module just
        /// read into this batch. May be less than `batch_size` near a file
        /// boundary, an event cap, or a shrunk final batch; 0 means unset.
        size_t events_in_batch{0};

        /// Raw pixel hits, keyed by detector name.
        std::map<std::string, std::vector<std::shared_ptr<Pixel>>> pixels;
        /// Raw digitizer waveforms, keyed by detector name.
        std::map<std::string, std::vector<std::shared_ptr<Waveform>>> waveforms;
        /// Clusters built from #pixels, keyed by detector name.
        std::map<std::string, std::vector<std::shared_ptr<Cluster>>> clusters;
        /// Tracks fitted from #clusters.
        std::vector<std::shared_ptr<Track>> tracks;
    };

    /**
     * @brief Base class for every reconstruction pipeline module.
     *
     * A module is constructed once from its own TOML config block, then
     * has initialize() called once, run() called once per batch (see
     * DataBatch), and finalize() called once after the last batch.
     */
    class Module {
    public:
        /**
         * @param cfg This module's own configuration block.
         * @param global_cfg The framework-wide `[Warlock]` configuration.
         * @param pool The thread pool shared by the whole pipeline.
         */
        Module(Configuration cfg, Configuration global_cfg, ThreadPool* pool);
        virtual ~Module() = default;

        /// Called once, before the first batch, in pipeline order.
        virtual void initialize() {}
        /// Called once per batch, in pipeline order (unless runsConcurrently()).
        virtual void run(DataBatch& batch) = 0;
        /// Called once, after the last batch, in pipeline order.
        virtual void finalize() {}

        /// Fraction of this module's own work completed so far, in [0,1],
        /// or a negative value if the module doesn't track progress itself
        /// (in which case the framework falls back to event-count progress).
        virtual double getProgress() const { return -1.0; }

        /**
         * @brief Whether this module may run concurrently with the rest of
         * the pipeline instead of blocking it.
         *
         * Reads the `run_concurrently` key from this module's own config
         * block, defaulting to `false` - meaning this module's run() fully
         * completes before the next module's run() starts, required for any
         * module that reads a DataBatch field another module writes.
         *
         * This is a config key rather than a per-class constant because
         * safety depends on which other modules share the pipeline, not on
         * this module's type alone: overlapping is safe only if no other
         * module in the same config reads or writes the DataBatch fields
         * this one touches (e.g. WaveformProcessingCAEN racing with
         * WaveformSelector over DataBatch::waveforms). Set
         * `run_concurrently = true` only after verifying, for that
         * specific config's module list, that nothing downstream touches
         * the same fields. FrameworkManager then launches the module
         * asynchronously and continues the pipeline immediately, so its
         * cost overlaps the rest of the batch instead of adding to it.
         */
        virtual bool runsConcurrently() const { return config.get<bool>("run_concurrently", false); }

        /// This module's configured instance name (e.g. `"WaveformSelector_CAEN_UZH_0"`).
        std::string getName() const;

        /// Sets this module's own output subdirectory (`<output_directory>/<name>/`).
        void setOutputPath(std::filesystem::path path);
        /// @return This module's own output subdirectory.
        std::filesystem::path getOutputPath() const;

    protected:
        Configuration config;        ///< This module's own configuration block.
        Configuration global_config; ///< The framework-wide `[Warlock]` configuration.
        ThreadPool* thread_pool;     ///< The thread pool shared by the whole pipeline.
        std::filesystem::path output_path_; ///< This module's own output subdirectory.
    };
}
#endif
