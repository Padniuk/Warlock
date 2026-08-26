/**
 * @file Logger.hpp
 * @brief Thread-safe logging utility for the Warlock framework.
 */

#ifndef WARLOCK_LOGGER_HPP
#define WARLOCK_LOGGER_HPP

#include <iostream>
#include <string>
#include <chrono>
#include <iomanip>
#include <mutex>
#include <functional>
#include <cstdio>
#ifndef _WIN32
#include <unistd.h>
#endif

namespace framework {
    /**
     * @enum LogLevel
     * @brief Defines the severity levels for log messages.
     */
    enum class LogLevel { TRACE, DEBUG, INFO, WARNING, ERROR, STATUS };

    /**
     * @brief Singleton, thread-safe, level-filtered logger with aligned,
     * colorized `|time| (LEVEL) [module] message` output.
     *
     * Use the WR_LOG macro to log from within a Module (it supplies the
     * module name via `getName()`); call log() directly anywhere else.
     */
    class Logger {
    public:
        /// Sets the minimum level that will actually print.
        static void set_level(LogLevel level) { get_instance().min_level_ = level; }

        /**
         * @brief Registers a callback run just before every message is
         * printed, while the logger's own mutex is held.
         *
         * Lets FrameworkManager erase its in-place progress bar line first,
         * so a log message can never land glued onto the end of the bar's
         * `\r`-updating line. Pass an empty std::function to clear it once
         * no bar is active.
         */
        static void setPreLogHook(std::function<void()> hook) {
            get_instance().pre_log_hook_ = std::move(hook);
        }

        /**
         * @brief Sets the column width every log line pads its `[module]`
         * field out to, so message text always starts at the same column
         * regardless of which module logged it. 0 (the default) disables
         * padding.
         */
        static void setModuleColumnWidth(size_t width) {
            get_instance().module_width_ = width;
        }

        /**
         * @brief Logs one message. Prefer the WR_LOG macro from within a Module.
         * @param level Severity/category of the message.
         * @param module Name shown in the `[module]` field.
         * @param message The message text.
         */
        static void log(LogLevel level, const std::string& module, const std::string& message) {
            auto& logger = get_instance();
            if (level < logger.min_level_) return;

            std::lock_guard<std::mutex> lock(logger.mutex_);
            if (logger.pre_log_hook_) logger.pre_log_hook_();

            auto now = std::chrono::system_clock::now();
            auto now_time_t = std::chrono::system_clock::to_time_t(now);
            auto now_ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

            // Brackets stay snug around the level word/module name itself
            // ("[STATUS]"/"{Framework}", not "[STATUS ]"/"{Framework   }")
            // - the padding needed to keep everything after them
            // column-aligned goes OUTSIDE the closing bracket/brace instead,
            // since padding baked inside is what made shorter words/names
            // look like they had a stray gap inside the tag itself.
            std::string word = to_string(level);
            size_t level_pad = (word.size() < kLevelWidth) ? (kLevelWidth - word.size()) : 0;
            size_t module_pad = (module.size() < logger.module_width_) ? (logger.module_width_ - module.size()) : 0;

            bool color = logger.colorEnabled();
            std::cout << "|" << std::put_time(std::localtime(&now_time_t), "%H:%M:%S")
                      << "." << std::setfill('0') << std::setw(3) << now_ms.count() << "| (";
            if (color) std::cout << levelColor(level);
            std::cout << word;
            if (color) std::cout << "\033[0m";
            std::cout << ")" << std::string(level_pad, ' ') << " [";
            if (color) std::cout << "\033[1m";
            std::cout << module;
            if (color) std::cout << "\033[0m";
            std::cout << "]" << std::string(module_pad, ' ') << " " << message << std::endl;
        }

    private:
        Logger() : min_level_(LogLevel::INFO) {}
        static Logger& get_instance() {
            static Logger instance;
            return instance;
        }

        // Cached once - colors/bold are actively unwanted noise once stdout
        // is redirected to a file or piped (log files with raw ANSI escapes
        // in them are a pain to grep/read later), so only emit them when
        // actually writing to an interactive terminal.
        static bool colorEnabled() {
#ifndef _WIN32
            static const bool tty = ::isatty(fileno(stdout)) != 0;
            return tty;
#else
            return false;
#endif
        }

        // "WARNING" is the longest level name - everything after the
        // closing "]" gets padded out to account for it, so that column
        // lands at the same place regardless of which level fired.
        static constexpr size_t kLevelWidth = 7;

        static std::string to_string(LogLevel level) {
            switch(level) {
                case LogLevel::TRACE:   return "TRACE";
                case LogLevel::DEBUG:   return "DEBUG";
                case LogLevel::INFO:    return "INFO";
                case LogLevel::WARNING: return "WARNING";
                case LogLevel::ERROR:   return "ERROR";
                case LogLevel::STATUS:  return "STATUS";
                default:                return "?";
            }
        }

        // Matches corry's convention of coloring by severity - dim for the
        // least important levels, up through bold red for ERROR; STATUS
        // (major run-level milestones, not a severity in the usual sense)
        // gets bold green so it reads as "notable/good" rather than urgent.
        static const char* levelColor(LogLevel level) {
            switch(level) {
                case LogLevel::TRACE:   return "\033[90m";   // bright black / gray
                case LogLevel::DEBUG:   return "\033[36m";   // cyan
                case LogLevel::INFO:    return "\033[37m";   // white
                case LogLevel::WARNING: return "\033[33m";   // yellow
                case LogLevel::ERROR:   return "\033[1;31m"; // bold red
                case LogLevel::STATUS:  return "\033[1;32m"; // bold green
                default:                return "\033[0m";
            }
        }

        LogLevel min_level_;
        std::mutex mutex_;
        std::function<void()> pre_log_hook_;
        size_t module_width_{0};
    };
}

/**
 * @brief Standard logging macro.
 * @param lvl The LogLevel (e.g., INFO, DEBUG, STATUS).
 * @param msg The message string.
 */
#define WR_LOG(lvl, msg) framework::Logger::log(framework::LogLevel::lvl, getName(), msg)

#endif
