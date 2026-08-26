/**
 * @file FrameworkManager.hpp
 * @brief Owns the reconstruction pipeline, the shared thread pool, and the
 * batch-processing run loop.
 */

#ifndef WARLOCK_FRAMEWORKMANAGER_HPP
#define WARLOCK_FRAMEWORKMANAGER_HPP

#include "core/Module.hpp"
#include "core/utils/ThreadPool.hpp"
#include "core/config/Configuration.hpp"
#include <vector>
#include <memory>
#include <string>

namespace framework {

    /**
     * @brief Top-level driver: parses a run's TOML config into a Module
     * pipeline, then repeatedly loads and processes batches of events
     * until the configured event or track limit is reached.
     */
    class FrameworkManager {
    public:
        FrameworkManager();

        /**
         * @brief Parses `config_file` into the module pipeline and the
         * framework-wide configuration.
         * @param config_file Path to the run's TOML configuration.
         * @param overrides Command-line `Module.key=value` overrides,
         * applied on top of the file before it's parsed into modules.
         */
        void loadConfiguration(const std::string& config_file, const std::vector<std::string>& overrides = {});

        /// Runs every module's initialize(), then the batch processing loop
        /// until completion, then every module's finalize().
        void run();

        std::string getName() const { return "Framework"; }

    private:
        /// Renders the single-line, in-place-updating progress bar.
        void printDynamicProgressBar(double progress, size_t events, size_t pixels, size_t clusters, size_t tracks, size_t waveforms);

        std::vector<std::unique_ptr<Module>> pipeline_; ///< Modules in declared (= scheduling) order.
        std::unique_ptr<ThreadPool> thread_pool_;        ///< Shared by every module in #pipeline_.
        Configuration global_config_;                    ///< The `[Warlock]` configuration block.

        /// Length of the last-printed progress line, so the next call can
        /// pad-blank any leftover tail with plain spaces instead of relying
        /// on the ANSI "erase entire line" (`\033[2K`) escape, which isn't
        /// reliably supported by every terminal - plain space-padding after
        /// a `\r` is about as universally supported as it gets.
        size_t prev_bar_len_{0};
    };
}
#endif
