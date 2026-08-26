/**
 * @file RawConverter.hpp
 * @brief Logic for converting EUDAQ2 raw data to Warlock HDF5 format.
 */

#ifndef WARLOCK_RAWCONVERTER_HPP
#define WARLOCK_RAWCONVERTER_HPP

#include <string>
#include <map>
#include <vector>
#include <chrono>

namespace framework {
    /** @brief Configuration for the raw converter tool. */
    struct ConverterConfig {
        std::string input_file;
        std::string output_file;
        size_t max_events = 0;
        // Raw EUDAQ2 producer descriptions to skip decoding entirely (e.g.
        // "CAEN_UZH" for a run where only telescope-plane hits are needed).
        // Matches corry's own discard_raw_events mechanism: the check runs
        // on the raw, undecoded (sub-)event, before the expensive
        // StdEventConverter::Convert() call, so a skipped CAEN digitizer's
        // per-sample correction tables (by far the dominant per-hit cost)
        // are never applied at all rather than applied and discarded.
        std::vector<std::string> skip_descriptions;
        // Raw EUDAQ2 producer descriptions that must be decoded strictly
        // one-at-a-time, PERMANENTLY, regardless of `threads` - for
        // converters whose thread-safety problem is an ongoing, per-event
        // one (e.g. CMSIT's converter compares each event's trigger ID
        // against the previous event's, so it needs strict sequential
        // order for its entire lifetime, not just at startup - see run()).
        std::vector<std::string> serial_descriptions;
        // Raw EUDAQ2 producer descriptions that only need sequential
        // decode for a one-time startup cost (lazy correction-table
        // loading, BORE-triggered Initialize()) - the first decode batch of
        // the file is forced fully sequential for these (see run()), so the
        // converter's Initialize()/table-load logic completes safely
        // exactly once; every batch after that, they're treated as safe to
        // parallelize like any other. The CAEN DT5742 converter is exactly
        // this case: its correction tables and per-device metadata are
        // read-only once loaded, and the per-sample correction math itself
        // is reentrant (X742CorrectionRoutines.c / X742DecodeRoutines.c use
        // only local variables and explicitly-passed pointers) - only the
        // lazy "have I loaded this yet" bookkeeping needs protecting.
        std::vector<std::string> warmup_descriptions;
        // Raw EUDAQ2 producer description that acts as the sole "clock" for
        // window formation - matches corry's own sync_by_event behavior
        // (EventLoaderEUDAQ2::is_within_event()), where whichever detector
        // loader is declared first in the .conf file (typically the
        // telescope's reference plane) defines each corryvreckan Event's
        // boundary from its own event number alone, and every other
        // producer's (sub-)event is matched against that boundary rather
        // than all producers being treated symmetrically. Left empty, the
        // converter falls back to symmetric "wait for all streams' numbers
        // to agree" windowing, which is wrong whenever the master
        // producer's number sequence isn't lockstep with everyone else's
        // (e.g. a rolling-shutter sensor whose frame period doesn't divide
        // evenly into the trigger period): a stray sub-event from any other
        // producer can spuriously open/close a window that corry never
        // forms, and the resulting offset compounds over the file. Set via
        // -m to whatever raw description the file uses for its anchor
        // producer (the converter logs every distinct description it sees
        // on first run).
        std::string master_description;
        // Worker threads used to parallelize StdEventConverter::Convert()
        // for descriptions not listed in serial_descriptions - defaults to
        // 1 (fully sequential decode) so a plain run with no flags is safe
        // for any producer mix; raise via -t once serial_descriptions
        // correctly lists every producer whose converter isn't verified
        // thread-safe.
        size_t threads = 1;
    };

    /**
     * @class RawConverter
     * @brief High-speed converter utilizing bulk-loading and sub-event filtering.
     */
    class RawConverter {
    public:
        explicit RawConverter(ConverterConfig cfg);
        void run();

    private:
        void printProgress(size_t current, size_t total, double elapsed_seconds, const std::string& diag_line = "");
        /// Blanks the current progress-bar line (via \r + spaces + \r,
        /// leaving the cursor at column 0 of a now-empty line) - call
        /// before any std::cout print that isn't itself a printProgress()
        /// call, so a mid-run diagnostic message can't collide with or
        /// leave behind a stale bar fragment. Same pattern as Warlock's own
        /// FrameworkManager::run() pre-log hook (prev_bar_len_).
        void eraseProgressLine();

        ConverterConfig cfg_;
        std::map<std::string, uint32_t> sensor_dict_;
        uint32_t next_sensor_id_ = 0;
        std::chrono::steady_clock::time_point start_time_;
        size_t prev_progress_len_ = 0; ///< Length of the last line printProgress() wrote - see eraseProgressLine().
    };
}

#endif