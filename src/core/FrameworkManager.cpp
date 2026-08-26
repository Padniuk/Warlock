#include "core/FrameworkManager.hpp"
#include "core/ModuleFactory.hpp"
#include "core/utils/PlotManager.hpp"
#include "core/utils/Logger.hpp"
#include "core/detector/GeometryManager.hpp"
#include <iostream>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <csignal>
#include <chrono>
#include <map>
#include <future>
#include <mutex>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/ioctl.h>
#include <unistd.h>
#endif

namespace {
    // Matches corry's compact progress-line style ("Ev: 107.8k Px: 6.28M
    // Tr: 282.9k") for large counters - 1 decimal at the thousands scale,
    // 2 decimals at the millions scale, plain otherwise.
    std::string formatCount(size_t value) {
        std::ostringstream oss;
        if (value >= 1000000) {
            oss << std::fixed << std::setprecision(2) << (static_cast<double>(value) / 1000000.0) << "M";
        } else if (value >= 1000) {
            oss << std::fixed << std::setprecision(1) << (static_cast<double>(value) / 1000.0) << "k";
        } else {
            oss << value;
        }
        return oss.str();
    }

    // Sizing the bar to the real terminal width avoids `\r` line-wrap
    // artifacts: if the line wraps, `\r` only returns to the start of the
    // wrapped row, not the true line start, leaving stale characters from a
    // previous, longer line. Falls back to a fixed width when there's no
    // real TTY to query (piped output, redirected to a file, etc.).
    int terminalWidth() {
        constexpr int fallback = 100;
#ifdef _WIN32
        CONSOLE_SCREEN_BUFFER_INFO csbi;
        if (GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi)) {
            return csbi.srWindow.Right - csbi.srWindow.Left + 1;
        }
#else
        struct winsize w;
        if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) == 0 && w.ws_col > 0) {
            return w.ws_col;
        }
#endif
        return fallback;
    }

    // `\033[?25l`/`\033[?25h` (hide/show cursor) are far more universally
    // supported than `\033[2K`/`\033[1A` - practically every terminal
    // implements cursor visibility. Without this, the block cursor sits
    // visibly at the end of the progress line for the whole run, since the
    // line is deliberately written without a trailing newline so the next
    // update can overwrite it in place.
    void hideCursor() { std::cout << "\033[?25l" << std::flush; }
    void showCursor() { std::cout << "\033[?25h" << std::flush; }

    // RAII: guarantees the cursor comes back on any ordinary exit path from
    // FrameworkManager::run() - normal completion, an early return, or a
    // C++ exception unwinding the stack. Does NOT cover a raw signal
    // termination (SIGINT/Ctrl+C default disposition kills the process
    // without unwinding, so no destructor runs) - see restoreCursorOnSigint
    // below for that case.
    struct CursorGuard {
        CursorGuard() { hideCursor(); }
        ~CursorGuard() { showCursor(); }
    };

#ifndef _WIN32
    // Only async-signal-safe calls belong in a signal handler: no iostream,
    // no std::string, no malloc, hence the raw write() to stdout.
    //
    // Calls _exit() rather than resetting to SIG_DFL and re-raising, since
    // POSIX delivers the signal to an arbitrary thread of the process (not
    // necessarily the one running the thread pool's owner), so that pattern
    // isn't reliable here. This loses the conventional "killed by signal"
    // wait-status a shell would otherwise see, so 128+signal (the standard
    // shell convention for that case) is used as the exit code instead.
    extern "C" void restoreCursorOnSigint(int sig) {
        const char* seq = "\033[?25h";
        ssize_t written = write(STDOUT_FILENO, seq, 6);
        (void)written; // nothing safe to do with a failed write here
        _exit(128 + sig);
    }
#endif
}

namespace framework {
    FrameworkManager::FrameworkManager() = default;

    void FrameworkManager::loadConfiguration(const std::string& config_file, const std::vector<std::string>& overrides) {
        toml::table tbl;
        try { tbl = toml::parse_file(config_file); }
        catch (const toml::parse_error& err) { throw std::runtime_error("TOML Parse Error: " + std::string(err.description())); }

        for (const auto& o : overrides) {
            auto dot_pos = o.find('.');
            auto eq_pos = o.find('=');
            if (dot_pos != std::string::npos && eq_pos != std::string::npos && eq_pos > dot_pos) {
                std::string mod = o.substr(0, dot_pos);
                std::string key = o.substr(dot_pos + 1, eq_pos - dot_pos - 1);
                std::string val = o.substr(eq_pos + 1);

                if (!tbl.contains(mod)) tbl.insert(mod, toml::table{});
                if (auto* mod_tbl = tbl.get_as<toml::table>(mod)) {
                    // Try the override value as genuine TOML first - lets a
                    // single -o flag set an array/number/bool natively (e.g.
                    // -o Reader.filenames=["a.h5","b.h5"]), not just a bare
                    // string. An unquoted bare value like a filename or path
                    // (run000268.h5, geometry/b7_geo1.toml) isn't valid TOML
                    // on its own, so this falls through to the plain
                    // string/bool handling below for exactly those cases.
                    bool applied = false;
                    try {
                        toml::table wrapper = toml::parse("v = " + val);
                        if (auto* node = wrapper.get("v")) {
                            mod_tbl->insert_or_assign(key, *node);
                            applied = true;
                        }
                    } catch (const toml::parse_error&) {
                        // Not valid standalone TOML - treat as a bare string/bool below.
                    }
                    if (!applied) {
                        if (val == "true") mod_tbl->insert_or_assign(key, true);
                        else if (val == "false") mod_tbl->insert_or_assign(key, false);
                        else mod_tbl->insert_or_assign(key, val);
                    }
                    std::cout << "[Warlock Config] OVERRIDE: " << mod << "." << key << " = " << val << "\n";
                }
            } else {
                std::cerr << "[Warlock Config] Warning: invalid override syntax. Use -o Module.key=val (" << o << ")\n";
            }
        }

        if (auto global = tbl.get_as<toml::table>("Warlock")) global_config_ = Configuration("Warlock", *global);
        else throw std::runtime_error("Configuration file must contain a [Warlock] block.");

        // Entries are parsed and the module-name log column width derived
        // from them before anything logs a line, so every line (including
        // the very first one) is already column-aligned.
        struct ModuleEntry {
            int line;
            std::string type_key;
            std::string instance_name;
            toml::table config_table;
        };
        std::vector<ModuleEntry> entries;
        for (auto&& [key, value] : tbl) {
            if (key == "Warlock") continue;
            std::string type_key = std::string(key.str());

            if (value.is_table()) {
                entries.push_back({static_cast<int>(value.source().begin.line), type_key, type_key, *value.as_table()});
            } else if (value.is_array()) {
                // [[ModuleType]] array-of-tables syntax: lets the same
                // module type be instantiated multiple times with different
                // configs (e.g. one WaveformSelector per detector). Each
                // element should carry its own "name" field so the instance
                // gets a distinguishable identity (histogram prefix, output
                // path, log lines); falls back to a bare index if absent.
                auto& arr = *value.as_array();
                int idx = 0;
                for (auto&& el : arr) {
                    if (!el.is_table()) { ++idx; continue; }
                    auto& sub = *el.as_table();
                    std::string instance_name = type_key + "_" + std::to_string(idx);
                    if (auto name_node = sub.get("name")) {
                        if (auto name_val = name_node->value<std::string>()) instance_name = type_key + "_" + *name_val;
                    }
                    entries.push_back({static_cast<int>(el.source().begin.line), type_key, instance_name, sub});
                    ++idx;
                }
            }
        }
        std::sort(entries.begin(), entries.end(),
                  [](const ModuleEntry& a, const ModuleEntry& b) { return a.line < b.line; });

        size_t module_col_width = getName().size();
        for (auto& e : entries) module_col_width = std::max(module_col_width, e.instance_name.size());
        Logger::setModuleColumnWidth(module_col_width);

        std::string level = global_config_.get<std::string>("log_level", "INFO");
        if(level == "DEBUG") Logger::set_level(LogLevel::DEBUG);
        else if(level == "TRACE") Logger::set_level(LogLevel::TRACE);

        std::string geo_path = global_config_.get<std::string>("detectors_file", "");
        if (!geo_path.empty()) {
            WR_LOG(INFO, "Loading geometry from: " + geo_path);
            GeometryManager::getInstance().loadGeometry(geo_path);
        }

        // "threads" is the key name every config in this project uses.
        size_t n_workers = global_config_.get<size_t>("threads", 4);
        thread_pool_ = std::make_unique<ThreadPool>(n_workers);

        WR_LOG(STATUS, "Building Reconstruction Pipeline (Source Order):");
        for (auto& e : entries) {
            WR_LOG(INFO, "  + Adding Module: " + e.instance_name);
            Configuration mod_cfg(e.instance_name, e.config_table);

            std::string base_out = global_config_.get<std::string>("output_directory", "output/");
            std::filesystem::path mod_dir = std::filesystem::path(base_out) / e.instance_name;

            auto mod = ModuleFactory::getInstance().createModule(e.type_key, mod_cfg, global_config_, thread_pool_.get());
            mod->setOutputPath(mod_dir);
            pipeline_.push_back(std::move(mod));
        }
    }

    void FrameworkManager::run() {
        // Hides the cursor for the duration of the dynamic progress bar
        // (restored automatically on any normal/exception exit via
        // CursorGuard's destructor), and makes sure Ctrl+C mid-run doesn't
        // leave it hidden afterward too - see CursorGuard and
        // restoreCursorOnSigint above for why both are needed.
        CursorGuard cursor_guard;
#ifndef _WIN32
        std::signal(SIGINT, restoreCursorOnSigint);
#endif

        // Erases the progress bar's current line before any log message
        // prints, so a WR_LOG call mid-run (e.g. an unrecognized-sensor
        // warning) can't land glued onto the end of the bar's
        // `\r`-updating line. Cleared at the end of this function once no
        // bar is active anymore.
        Logger::setPreLogHook([this]() {
            if (prev_bar_len_ > 0) {
                std::cout << "\r" << std::string(prev_bar_len_, ' ') << "\r";
                prev_bar_len_ = 0;
            }
        });

        // Real wall-clock elapsed time for the whole run (initialize()
        // through the final plot save), printed alongside the per-module
        // breakdown below. Not the same as summing that breakdown: once a
        // module runs concurrently (see Module::runsConcurrently()), its
        // timer window overlaps whatever runs after it, so the per-module
        // sum can exceed this real total.
        auto pipeline_start = std::chrono::steady_clock::now();

        for (auto& mod : pipeline_) mod->initialize();

        auto events_arr = global_config_.getArray<size_t>("number_of_events");
        size_t event_limit = 0;
        for (auto ev : events_arr) event_limit += ev;

        // Matches corry's global `number_of_tracks`: a single
        // framework-wide track-count target, independent of which specific
        // module produces or consumes tracks. 0 means unlimited.
        size_t track_limit = global_config_.get<size_t>("number_of_tracks", 0);
        bool has_track_target = track_limit > 0;

        size_t b_size = global_config_.get<size_t>("batch_size", 5000);
        size_t events_processed = 0;
        size_t pixels_processed = 0;
        size_t clusters_processed = 0;
        size_t tracks_processed = 0;
        size_t waveforms_processed = 0;

        // Per-module wall-clock time, accumulated across every batch -
        // matches corry's own end-of-run "Wall-clock timing" breakdown,
        // printed below.
        std::map<std::string, double> module_time_s;
        std::mutex module_time_mtx; // written from concurrent modules' own threads too

        // Printed once, before the first batch even starts - without this,
        // the bar's first appearance is only after batch 1 has been fully
        // processed, which for a slow first batch reads as the program
        // having hung at startup.
        printDynamicProgressBar(0.0, 0, 0, 0, 0, 0);

        while (true) {
            DataBatch batch;
            batch.batch_id = events_processed;

            // Shrink the next batch as a `number_of_tracks` target gets
            // close, so the run doesn't cluster/track a full extra chunk of
            // events beyond what's actually needed - that extra work is
            // never done by corry's per-event loop, and it would otherwise
            // overshoot the requested track count and skew an
            // apples-to-apples comparison against corry's output. The
            // event-limit-only case doesn't need this: EventLoaderHDF5
            // already caps precisely per file via `number_of_events`.
            if (has_track_target && events_processed > 0 && tracks_processed > 0 && tracks_processed < track_limit) {
                double tracks_per_event = static_cast<double>(tracks_processed) / static_cast<double>(events_processed);
                if (tracks_per_event > 0.0) {
                    size_t remaining_tracks = track_limit - tracks_processed;
                    // Undershoot the naive estimate: the running average
                    // rate won't exactly match this batch's actual yield,
                    // so requesting only 90% keeps an above-average batch
                    // from overshooting the target, at the cost of one more
                    // (much smaller) pass through the loop.
                    constexpr double safety_factor = 0.9;
                    size_t estimated_events = static_cast<size_t>(std::ceil(remaining_tracks / tracks_per_event * safety_factor));
                    batch.requested_batch_size = std::max<size_t>(estimated_events, 1);

                    // Once the remaining gap is small, drop to single-event
                    // batches - the finest granularity available - so the
                    // final overshoot is bounded by one event's yield,
                    // matching corry's own per-event stopping check. The
                    // margin is scaled by tracks_per_event rather than a
                    // fixed event count, since a single batch's local yield
                    // can run well above the cumulative average.
                    constexpr double single_event_margin_events = 20.0;
                    if (static_cast<double>(remaining_tracks) <= tracks_per_event * single_event_margin_events) {
                        batch.requested_batch_size = 1;
                    }
                }
            }

            // Modules that opted into runsConcurrently() (see
            // Module::runsConcurrently()) are launched here and left
            // running while the loop continues to the sequential modules
            // after them, instead of blocking on them. Their own chunked
            // work still goes through the same shared thread_pool_, so this
            // interleaves work rather than carving out dedicated cores.
            // Joined once the sequential pass finishes.
            std::vector<std::future<void>> concurrent_futures;
            for (auto& mod : pipeline_) {
                if (mod->runsConcurrently()) {
                    concurrent_futures.push_back(std::async(std::launch::async, [&mod, &batch, &module_time_s, &module_time_mtx]() {
                        auto t0 = std::chrono::steady_clock::now();
                        mod->run(batch);
                        double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                        std::lock_guard<std::mutex> lock(module_time_mtx);
                        module_time_s[mod->getName()] += dt;
                    }));
                    continue;
                }
                auto t0 = std::chrono::steady_clock::now();
                mod->run(batch);
                module_time_s[mod->getName()] +=
                    std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                if (batch.halt_pipeline) break;
            }
            // Must finish before `batch` goes out of scope at the end of
            // this loop iteration - the lambdas above captured it by
            // reference.
            for (auto& f : concurrent_futures) f.get();

            bool end_of_data = batch.pixels.empty() && batch.waveforms.empty() && batch.tracks.empty() &&
                               batch.clusters.empty();

            for (auto const& entry : batch.pixels) pixels_processed += entry.second.size();
            for (auto const& entry : batch.clusters) clusters_processed += entry.second.size();
            tracks_processed += batch.tracks.size();
            // Waveforms are independent of the events -> pixels -> clusters
            // -> tracks chain above (WaveformProcessingCAEN works straight off
            // each event's raw CAEN samples) - counted the same way regardless.
            for (auto const& entry : batch.waveforms) waveforms_processed += entry.second.size();
            // Counted only once we know this wasn't just the empty sentinel
            // batch signalling end-of-data, and using the event source's
            // own reported count (not the configured batch_size) since a
            // shrunk final batch reads fewer events than that.
            if (!end_of_data) events_processed += (batch.events_in_batch > 0 ? batch.events_in_batch : b_size);

            // `number_of_tracks`, when set, drives both the displayed
            // progress and the stopping decision in preference to the event
            // count - matching corry, where a track target always takes
            // priority over the event count. Modules that need to stop
            // consuming tracks precisely mid-batch (e.g. alignment modules)
            // still do so themselves via batch.halt_pipeline; this
            // framework-level check is the generic fallback for any
            // pipeline, with or without such a module.
            if (batch.halt_pipeline || end_of_data) {
                printDynamicProgressBar(1.0, events_processed, pixels_processed, clusters_processed, tracks_processed, waveforms_processed);
                std::cout << "\n\n";
                break;
            }

            double progress;
            if (has_track_target) {
                progress = static_cast<double>(tracks_processed) / track_limit;
            } else if (event_limit > 0) {
                progress = static_cast<double>(events_processed) / event_limit;
            } else {
                // Neither a track target nor an event count was given - fall
                // back to whatever generic progress a module can report
                // (e.g. an event loader's fraction-of-file-read).
                progress = 0.0;
                for (auto& mod : pipeline_) progress = std::max(progress, mod->getProgress());
            }

            printDynamicProgressBar(progress, events_processed, pixels_processed, clusters_processed, tracks_processed, waveforms_processed);

            if (has_track_target && tracks_processed >= track_limit) {
                printDynamicProgressBar(1.0, events_processed, pixels_processed, clusters_processed, tracks_processed, waveforms_processed);
                std::cout << "\n\n";
                break;
            }
            if (!has_track_target && event_limit > 0 && events_processed >= event_limit) {
                printDynamicProgressBar(1.0, events_processed, pixels_processed, clusters_processed, tracks_processed, waveforms_processed);
                std::cout << "\n\n";
                break;
            }
        }

        for (auto& mod : pipeline_) mod->finalize();

        // Matches corry's own end-of-run "Wall-clock timing" section -
        // sorted slowest-first (not pipeline order), so the dominant cost
        // is obvious at a glance.
        {
            std::vector<std::pair<std::string, double>> sorted_times(module_time_s.begin(), module_time_s.end());
            std::sort(sorted_times.begin(), sorted_times.end(),
                      [](auto const& a, auto const& b) { return a.second > b.second; });

            std::ostringstream ts;
            ts << "Wall-clock timing (seconds):\n";
            for (auto const& [name, secs] : sorted_times) {
                double ms_per_evt = (events_processed > 0) ? (secs * 1000.0 / static_cast<double>(events_processed)) : 0.0;
                ts << "  " << std::left << std::setw(24) << name << std::right
                   << std::fixed << std::setprecision(3) << std::setw(10) << secs << "s = "
                   << std::setprecision(4) << ms_per_evt << "ms/evt\n";
            }
            double total_pipeline_s = std::chrono::duration<double>(std::chrono::steady_clock::now() - pipeline_start).count();
            ts << "  " << std::left << std::setw(24) << "TOTAL (wall-clock)" << std::right
               << std::fixed << std::setprecision(3) << std::setw(10) << total_pipeline_s << "s"
               << " - real elapsed time; NOT the sum of the rows above (those can overlap under concurrency)\n";
            WR_LOG(STATUS, ts.str());
        }

        std::string plots_filename = global_config_.get<std::string>("output_plots_file", "plots.h5");
        std::string base_out = global_config_.get<std::string>("output_directory", "output/");
        std::filesystem::path full_plots_path = std::filesystem::path(base_out) / plots_filename;

        PlotManager::getInstance().saveToFile(full_plots_path.string());
        WR_LOG(STATUS, "Saved all module plots to: " + full_plots_path.string());

        Logger::setPreLogHook(nullptr);
    }

    void FrameworkManager::printDynamicProgressBar(double progress, size_t events, size_t pixels, size_t clusters, size_t tracks, size_t waveforms) {
        if (progress < 0.0) progress = 0.0;
        if (progress > 1.0) progress = 1.0;

        std::ostringstream prefix;
        prefix << std::setfill(' ') << std::setw(3) << static_cast<int>(progress * 100) << "% [";

        std::ostringstream suffix;
        // Ev -> Px -> Cl -> Tr is the actual production order (pixels come
        // from events, clusters from pixels, tracks from clusters). Wf is
        // listed last because it genuinely isn't part of that chain -
        // WaveformProcessingCAEN works straight off each event's raw CAEN
        // samples, independent of pixels/clustering/tracking.
        suffix << "]  Ev: " << formatCount(events) << "  Px: " << formatCount(pixels)
               << "  Cl: " << formatCount(clusters) << "  Tr: " << formatCount(tracks)
               << "  Wf: " << formatCount(waveforms);

        std::string prefix_s = prefix.str();
        std::string suffix_s = suffix.str();

        // Size the bar to whatever room is left after the percentage prefix
        // and stats suffix, so the whole line fits the actual terminal
        // width - clamped so it's never so narrow it stops looking like a
        // bar, and never wider than the original 40.
        constexpr int min_width = 10, max_width = 40;
        int available = terminalWidth() - static_cast<int>(prefix_s.size() + suffix_s.size());
        int width = std::clamp(available, min_width, max_width);
        int pos = static_cast<int>(width * progress);

        std::ostringstream line;
        line << prefix_s;
        for (int i = 0; i < width; ++i) line << (i < pos ? "█" : " ");
        line << suffix_s;

        // Single line, overwritten via a leading \r, padded with plain
        // spaces to blank out any leftover tail from a longer previous
        // line - deliberately NOT using `\033[2K` (erase line) or `\033[1A`
        // (cursor up), which aren't reliably supported by every terminal.
        //
        // display_len is computed from prefix+width+suffix rather than
        // s.size(): "█" is a 3-byte UTF-8 sequence but occupies 1 terminal
        // column, so byte length overcounts the real on-screen width and
        // can leave stale characters un-blanked when the line shrinks.
        std::string s = line.str();
        size_t display_len = prefix_s.size() + static_cast<size_t>(width) + suffix_s.size();
        std::cout << "\r" << s;
        if (display_len < prev_bar_len_) std::cout << std::string(prev_bar_len_ - display_len, ' ');
        std::cout << std::flush;
        prev_bar_len_ = display_len;
    }
}
