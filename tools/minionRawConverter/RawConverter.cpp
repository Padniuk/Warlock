/**
 * @file RawConverter.cpp
 * @brief RawConverter implementation - see RawConverter.hpp for the class/config docs.
 */

#include "RawConverter.hpp"
#include "AppendableDataset.hpp"
#include "core/utils/ThreadPool.hpp"

#include <eudaq/FileReader.hh>
#include <eudaq/StandardEvent.hh>
#include <eudaq/StdEventConverter.hh>

#include <iostream>
#include <iomanip>
#include <sstream>
#include <queue>
#include <optional>
#include <algorithm>
#include <thread>
#include <set>
#include <map>
#include <regex>
#include <cctype>
#include <cmath>

namespace framework {

    namespace {
        // One producer's contribution to a single physical trigger. A raw
        // EUDAQ2 event read from the file may bundle several of these as
        // sub-events (multi-board setups), or arrive unbundled as a single
        // StandardEvent - both cases get flattened into this same unit
        // before synchronization, so the rest of the converter doesn't need
        // to care which form the file happened to store them in.
        struct AtomicEvent {
            std::shared_ptr<const eudaq::StandardEvent> evt;
            uint32_t event_number;
            // True if this (sub-)event's raw description matched
            // ConverterConfig::master_description - see the windowing loop
            // in RawConverter::run() for how this drives window formation.
            bool is_master = false;
        };

        // CAEN DT5742 ADC<->Volts calibration constants and helpers, mirrored
        // from eudaq/user/caen-dt5742's CAENDT5742RawEvent2StdEventConverter.cc
        // (ADC12ToVolts / DACToBaselineVolts / ParseChannelDACMap - same
        // constants, same regex, same default DAC of 32768 for an
        // unconfigured channel). Used to invert that plugin's calibration:
        // eudaq's converter only exposes the fully-calibrated Volts waveform
        // via GetWaveform(), with no public accessor for the intermediate
        // DRS4-corrected ADC code. Since the final step is the linear
        // function `v = (a - 2048) * (Vpp/4096) + dc_offset`, inverting it
        // recovers that ADC code exactly (modulo float rounding), letting
        // `samples` store it as int16 (2 bytes) instead of Volts as
        // double/float (8/4 bytes) - see the `ds_wf_samples` comment below.
        constexpr double kCaenAdcMid = 2048.0;
        constexpr double kCaenAdcLsbDen = 4096.0;
        constexpr double kCaenVpp = 1.0;  // DT5742 default input dynamic range

        double caenDacToBaselineVolts(uint32_t dac) {
            return (static_cast<double>(dac) / 65535.0 - 0.5) * kCaenVpp;
        }

        std::map<int, uint32_t> parseCaenChannelDacMap(const std::string& text) {
            std::map<int, uint32_t> result;
            static const std::regex re(R"(CH(\d+)\s*:\s*(\d+))");
            for (std::sregex_iterator it(text.begin(), text.end(), re); it != std::sregex_iterator(); ++it) {
                const std::smatch m = *it;
                result[std::stoi(m[1].str())] = static_cast<uint32_t>(std::stoul(m[2].str()));
            }
            return result;
        }
    }

    RawConverter::RawConverter(ConverterConfig cfg) : cfg_(std::move(cfg)) {}

    void RawConverter::eraseProgressLine() {
        if (prev_progress_len_ > 0) {
            std::cout << "\r" << std::string(prev_progress_len_, ' ') << "\r";
            prev_progress_len_ = 0;
        }
    }

    void RawConverter::printProgress(size_t current, size_t total, double elapsed_seconds, const std::string& diag_line) {
        // Remembers the length of the last line printed (prev_progress_len_,
        // a class member so eraseProgressLine() can also use it - see its
        // docs) so this call can pad out to at least that length, blanking
        // any leftover tail from a longer previous line without relying on
        // an erase-to-end-of-line escape code. Everything (bar + diagnostics)
        // lives on one line, overwritten via a leading \r, since cursor-up
        // escape sequences don't render reliably in every terminal.
        std::ostringstream line;
        double ev_per_sec = (elapsed_seconds > 0) ? (current / elapsed_seconds) : 0;
        line << " [minion] ";
        int bar_chars = 0; ///< Count of "█" glyphs actually written below - see display_len's docs.
        if (total > 0) {
            int barWidth = 40;
            float progress = static_cast<float>(current) / static_cast<float>(total);
            int pos = static_cast<int>(barWidth * progress);
            line << std::setw(3) << static_cast<int>(progress * 100) << "% [";
            for (int i = 0; i < barWidth; ++i) line << (i < pos ? "█" : " ");
            line << "] " << current << "/" << total;
            bar_chars = pos;
        } else {
            line << current << " events";
        }
        line << " [" << std::fixed << std::setprecision(1) << ev_per_sec << " ev/s]";
        if (!diag_line.empty()) line << "  |  " << diag_line;

        // Tracks/compares DISPLAY WIDTH (terminal columns), not
        // std::string::size() (bytes): "█" is a 3-byte UTF-8 sequence but
        // only ever occupies 1 column, so byte length overcounts by 2
        // bytes per filled bar segment. Same fix as Warlock's own
        // FrameworkManager::printDynamicProgressBar().
        std::string s = line.str();
        size_t display_len = s.size() - static_cast<size_t>(2 * bar_chars);
        std::cout << "\r" << s;
        if (display_len < prev_progress_len_) std::cout << std::string(prev_progress_len_ - display_len, ' ');
        std::cout << std::flush;
        prev_progress_len_ = display_len;
    }

    void RawConverter::run() {
        auto reader = eudaq::Factory<eudaq::FileReader>::MakeUnique(eudaq::str2hash("native"), cfg_.input_file);

        // File-creation tuning: the default B-tree branching factor for
        // indexing chunks in an extensible dataset (istore_k = 32) is sized
        // for modest chunk counts. A multi-million-event run fills thousands
        // of chunks per dataset, so a wider branching factor keeps that
        // index shallow instead of letting per-write lookup cost grow with
        // file size.
        H5::FileCreatPropList fcpl;
        fcpl.setIstorek(128);

        // File-level raw-data chunk cache defaults, as a safety net for any
        // dataset opened without its own explicit DAPL (see
        // AppendableDataset.cpp for the per-dataset settings that matter
        // most). mdc_nelmts is a legacy parameter HDF5 ignores since the
        // metadata cache became self-tuning; -1 keeps the library default.
        H5::FileAccPropList facpl;
        facpl.setCache(-1, 10007, 64 * 1024 * 1024, 0.75);

        H5::H5File h5file(cfg_.output_file, H5F_ACC_TRUNC, fcpl, facpl);
        H5::Group grp_ev = h5file.createGroup("/events");
        H5::Group grp_hi = h5file.createGroup("/hits");

        using converter::AppendableDataset;
        AppendableDataset ds_ev_id(grp_ev, "event_id", H5::PredType::NATIVE_UINT64);
        AppendableDataset ds_ev_num(grp_ev, "event_number", H5::PredType::NATIVE_UINT32);
        AppendableDataset ds_hi_evid(grp_hi, "event_id", H5::PredType::NATIVE_UINT64);
        AppendableDataset ds_hi_col(grp_hi, "col", H5::PredType::NATIVE_INT32);
        AppendableDataset ds_hi_row(grp_hi, "row", H5::PredType::NATIVE_INT32);
        AppendableDataset ds_hi_sens(grp_hi, "sensor_id", H5::PredType::NATIVE_UINT32);
        AppendableDataset ds_hi_plane(grp_hi, "plane_id", H5::PredType::NATIVE_UINT32);
        AppendableDataset ds_hi_ts(grp_hi, "timestamp", H5::PredType::NATIVE_DOUBLE);
        AppendableDataset ds_hi_haswf(grp_hi, "has_waveform", H5::PredType::NATIVE_UINT8);
        // Raw hardware channel (0-17 for a CAEN DT5742 board), parsed from
        // the CAEN plugin's own per-pixel aux-info string
        // ("<dutname>:CH<N>:col<C>:row<R>", see extractChannel() below) -
        // -1 for any hit where that string isn't present/parseable (every
        // non-CAEN sensor). Lets downstream geometry split one physical
        // multi-die CAEN board's raw hits by die (or, in principle, by
        // trigger group) without needing separate raw sensor identities.
        AppendableDataset ds_hi_channel(grp_hi, "channel", H5::PredType::NATIVE_INT32);
        // Waveform datasets are compact: one row per hit that actually HAS
        // waveform data (see has_waveform), not one row per hit overall.
        // Most detectors in a mixed telescope+DUT run (e.g. MIMOSA26) never
        // produce waveforms, so padding every hit with 1024 zeroed doubles
        // would waste most of both conversion time and file size.
        //
        // Stored as int16, holding the DRS4-corrected ADC code (~12 bits,
        // centered on 2048 - see kCaenAdcMid above), not the calibrated
        // Volts value: DRS4 correction is expensive/stateful (per-device
        // correction-table files, not thread-safe on first load - see
        // EnsureCorrectionTablesLoaded) and must happen here, once, while
        // the final ADC->Volts step is a trivial stateless linear formula
        // deferred to WaveformProcessingCAEN (run at analysis time, only on
        // waveforms actually read). append2D() writes from an in-memory
        // double buffer regardless of this dataset's on-disk type (values
        // already rounded to integers, so the double->int16 write is an
        // exact truncation); HDF5 converts on both write and read, so reads
        // still come back as NATIVE_DOUBLE - only the physical meaning
        // (ADC code, not Volts) has changed. See WaveformProcessingCAEN.cpp's
        // run() for the corresponding decode.
        AppendableDataset ds_wf_samples(grp_hi, "samples", H5::PredType::NATIVE_INT16, 2);
        AppendableDataset ds_wf_t0(grp_hi, "wf_t0", H5::PredType::NATIVE_DOUBLE);
        AppendableDataset ds_wf_dt(grp_hi, "wf_dt", H5::PredType::NATIVE_DOUBLE);

        std::vector<uint64_t> b_ev_id, b_hi_evid;
        std::vector<uint32_t> b_ev_num, b_hi_sens, b_hi_plane;
        std::vector<int32_t> b_hi_col, b_hi_row, b_hi_channel;
        std::vector<double> b_hi_ts;
        std::vector<uint8_t> b_hi_haswf;
        std::vector<double> b_wf_samples, b_wf_t0, b_wf_dt;

        // Per-channel DC-offset baseline (Volts), needed to invert eudaq's
        // ADC12ToVolts and recover the corrected-ADC-code value stored in
        // `samples` (see kCaenAdcMid comment above). caen_dac_map is parsed
        // once from the CAEN BORE event's "channel_dc_offset_dac_map" tag,
        // if present; caen_channel_dc_offset_volts memoizes the resolved
        // Volts value per channel (falling back to eudaq's default DAC of
        // 32768, mid-scale/"no offset", for any channel missing from the
        // tag/map), and is persisted as `/hits` attributes at the end so
        // WaveformProcessingCAEN can reverse the same calibration without
        // re-parsing raw-file tags itself.
        std::map<int, uint32_t> caen_dac_map;
        std::map<int32_t, double> caen_channel_dc_offset_volts;
        auto caenDcOffsetVolts = [&](int32_t channel) -> double {
            auto it = caen_channel_dc_offset_volts.find(channel);
            if (it != caen_channel_dc_offset_volts.end()) return it->second;
            uint32_t dac = 32768u;
            auto dac_it = caen_dac_map.find(channel);
            if (dac_it != caen_dac_map.end()) dac = dac_it->second;
            double dc = caenDacToBaselineVolts(dac);
            caen_channel_dc_offset_volts[channel] = dc;
            return dc;
        };

        uint64_t global_id = 0;
        start_time_ = std::chrono::steady_clock::now();

        // Throughput diagnostics: per-window hit/waveform counts and a
        // decode-time-vs-wall-time split, printed alongside the normal
        // progress line every 100 events.
        uint64_t diag_hits = 0, diag_wf_hits = 0, diag_subevents_decoded = 0;
        double diag_decode_ms = 0.0;
        uint64_t diag_prev_global_id = 0, diag_prev_hits = 0, diag_prev_wf = 0, diag_prev_subevents = 0;
        double diag_prev_decode_ms = 0.0, diag_prev_elapsed_ms = 0.0;

        std::cout << "[minion] Converting"
                  << (cfg_.max_events > 0 ? (" up to " + std::to_string(cfg_.max_events) + " events") : "")
                  << "..." << std::endl;

        // Raw, UNDECODED (sub-)events read from the file but not yet
        // decoded/consumed - populated a whole raw eudaq::Event at a time
        // (which may unpack into several of these), drained one at a time.
        std::queue<eudaq::EventSPC> pending;
        bool reader_eof = false;

        // Reading the raw file is inherently sequential (one binary stream),
        // and the event-number windowing below needs (sub-)events in their
        // original order too - but StdEventConverter::Convert() itself is,
        // for a WELL-BEHAVED plugin, a pure function of one (sub-)event,
        // independent of any others, so it can in principle run across a
        // thread pool over a batch of raw (sub-)events, with the batch fed
        // through the windowing logic below afterward in original order.
        //
        // Not all plugins are well-behaved, though: the CAEN DT5742 converter
        // (CAENDT5742RawEvent2StdEventConverter.cc) keeps its correction
        // tables and per-device metadata in `static` class member maps
        // (_x742_correction_tables, _x742_correction_table_loaded, _name,
        // _dut_channel_arrangement, ...), and EnsureCorrectionTablesLoaded()
        // does an unsynchronized check-then-load-then-write on those maps
        // the first time each device/group is seen - decoding two CAEN
        // (sub-)events concurrently on their first appearance is a genuine
        // data race. eudaq/user/eudet's MIMOSA26 converter
        // (NiRawEvent2StdEventConverter.cc), by contrast, only has
        // compile-time `static const` members and is safe to parallelize.
        //
        // Since a single raw file mixes producers and this tool has no
        // reliable way to know every plugin's thread-safety in general,
        // decoding defaults to sequential (threads=1) - safe for any run,
        // any producer mix, with no flag required. Parallel decode is
        // opt-in only, via -t, for someone who has verified every producer
        // type in their specific file is safe to decode concurrently.
        size_t n_threads = cfg_.threads > 0 ? cfg_.threads : 1;
        ThreadPool pool(n_threads);
        // Sub-events (not combined events) read ahead before the batch gets
        // decoded. Kept fairly small rather than maximizing thread-dispatch
        // amortization, because warmup_descriptions (see below) only lets a
        // producer into the parallel pool from the second batch onward - a
        // large first batch means a large fully-sequential prefix before
        // that benefit kicks in.
        constexpr size_t kDecodeBatch = 4000;

        std::vector<eudaq::EventSPC> raw_batch;
        std::vector<std::optional<AtomicEvent>> decoded_batch;
        size_t decoded_idx = 0;

        // Logged once per distinct raw event description encountered, so a
        // first run with no -s/-x tells you exactly what strings this file
        // actually uses, instead of guessing.
        std::set<std::string> seen_descriptions;

        auto fill_raw_batch = [&]() {
            raw_batch.clear();
            while (raw_batch.size() < kDecodeBatch && !reader_eof) {
                while (pending.empty() && !reader_eof) {
                    auto raw = reader->GetNextEvent();
                    if (!raw) { reader_eof = true; break; }

                    auto subs = raw->GetSubEvents();
                    if (subs.empty()) pending.push(raw);
                    else for (const auto& sub : subs) pending.push(sub);
                }
                if (pending.empty()) break;

                auto raw_sub = pending.front();
                pending.pop();

                const std::string& desc = raw_sub->GetDescription();
                if (seen_descriptions.insert(desc).second) {
                    // Blank the active progress-bar line first (see
                    // eraseProgressLine()'s docs) so this print doesn't
                    // collide with or leave behind a stale bar fragment.
                    eraseProgressLine();
                    std::cout << "[minion] New raw event description: \"" << desc << "\"" << std::endl;
                }

                // Same tag eudaq's own CAEN plugin reads on its BORE event
                // (see caen_dac_map's declaration above) - read here too,
                // straight off the raw sub-event, since minionRawConverter
                // needs it independently to invert the plugin's calibration.
                if (desc == "CAENDT5742" && raw_sub->IsBORE() && caen_dac_map.empty()) {
                    caen_dac_map = parseCaenChannelDacMap(raw_sub->GetTag("channel_dc_offset_dac_map"));
                }

                // BORE/EORE (sub-)events are deliberately NOT filtered out
                // here, even though they carry no hit data themselves:
                // stateful converter plugins (e.g. CAEN's) only run their
                // Initialize() - which populates correction tables and
                // per-device metadata - as a side effect of Convert() being
                // called with a BORE event. Filtering BORE out before ever
                // calling Convert() would skip that initialization, and
                // every subsequent hit from that producer would fail to
                // decode. Convert() itself correctly rejects BORE/EORE as
                // containing no real hit data (returns false, or a decoded
                // event with zero planes), which decode_batch_parallel below
                // treats as "nothing to contribute", so letting them through
                // here is free of side effects beyond the Initialize() we
                // want.
                if (!cfg_.skip_descriptions.empty()) {
                    if (std::find(cfg_.skip_descriptions.begin(), cfg_.skip_descriptions.end(), desc) !=
                        cfg_.skip_descriptions.end()) {
                        continue;
                    }
                }

                raw_batch.push_back(std::move(raw_sub));
            }
        };

        // Set once the very first decode batch has been processed. Producers
        // in warmup_descriptions are only forced serial while this is still
        // false - BORE events (and hence Initialize()/lazy table-loading)
        // occur at the very start of the file, so the first batch is enough
        // to let a converter like CAEN's finish its one-time setup safely;
        // every batch after that treats it as any other producer.
        bool warmed_up = false;

        auto decode_batch_parallel = [&]() {
            decoded_batch.assign(raw_batch.size(), std::nullopt);
            if (raw_batch.empty()) return;

            auto decode_one = [&](size_t i) {
                auto decoded = eudaq::StandardEvent::MakeShared();
                if (eudaq::StdEventConverter::Convert(raw_batch[i], decoded, nullptr)) {
                    // Accept a decoded (sub-)event if it carries real pixel
                    // data OR a real StandardEvent-level timestamp. Of the
                    // four producers handled here, only TLU calls
                    // SetTimeBegin()/SetTimeEnd(); CAEN, NI (MIMOSA26), and
                    // CMSIT never do, and none push a nonzero per-pixel
                    // timestamp either. TLU itself never produces a
                    // StandardPlane (NumPlanes() is always 0), so requiring
                    // NumPlanes() > 0 here would silently drop TLU's decoded
                    // result on every event, discarding the only real trigger
                    // clock in the file - see close_window() below for how
                    // that timestamp is applied across the window.
                    bool has_planes = decoded->NumPlanes() > 0;
                    bool has_time = decoded->GetTimeBegin() != 0 || decoded->GetTimeEnd() != 0;
                    if (has_planes || has_time) {
                        bool is_master = !cfg_.master_description.empty() &&
                                          raw_batch[i]->GetDescription() == cfg_.master_description;
                        decoded_batch[i] = AtomicEvent{decoded, decoded->GetEventNumber(), is_master};
                    }
                }
            };

            // Route each (sub-)event by producer: serial_descriptions
            // decode one at a time forever; warmup_descriptions do too, but
            // only until the first batch has gone through. Everything else
            // is fair game for the pool.
            std::vector<size_t> parallel_idx;
            parallel_idx.reserve(raw_batch.size());
            for (size_t i = 0; i < raw_batch.size(); ++i) {
                bool force_serial = false;
                if (!cfg_.serial_descriptions.empty() || (!warmed_up && !cfg_.warmup_descriptions.empty())) {
                    const auto& desc = raw_batch[i]->GetDescription();
                    if (std::find(cfg_.serial_descriptions.begin(), cfg_.serial_descriptions.end(), desc) !=
                        cfg_.serial_descriptions.end()) {
                        force_serial = true;
                    } else if (!warmed_up &&
                               std::find(cfg_.warmup_descriptions.begin(), cfg_.warmup_descriptions.end(), desc) !=
                                   cfg_.warmup_descriptions.end()) {
                        force_serial = true;
                    }
                }
                if (force_serial) decode_one(i);
                else parallel_idx.push_back(i);
            }

            if (n_threads > 1 && !parallel_idx.empty()) {
                size_t chunk = (parallel_idx.size() + n_threads - 1) / n_threads;
                std::vector<std::future<void>> futures;
                for (size_t t = 0; t < n_threads; ++t) {
                    size_t start = t * chunk, end = std::min(start + chunk, parallel_idx.size());
                    if (start >= end) continue;
                    futures.push_back(pool.submit([&, start, end]() {
                        for (size_t k = start; k < end; ++k) decode_one(parallel_idx[k]);
                    }));
                }
                for (auto& f : futures) f.get();
            } else {
                for (size_t i : parallel_idx) decode_one(i);
            }

            warmed_up = true;
        };

        auto next_atomic = [&]() -> std::optional<AtomicEvent> {
            while (true) {
                if (decoded_idx >= decoded_batch.size()) {
                    fill_raw_batch();
                    if (raw_batch.empty()) return std::nullopt;
                    diag_subevents_decoded += raw_batch.size();
                    auto diag_t0 = std::chrono::steady_clock::now();
                    decode_batch_parallel();
                    diag_decode_ms += std::chrono::duration<double, std::milli>(
                        std::chrono::steady_clock::now() - diag_t0).count();
                    decoded_idx = 0;
                }
                auto& opt = decoded_batch[decoded_idx++];
                if (opt) return std::move(opt);
                // Decode failed for this (sub-)event - try the next one.
            }
        };

        // The CAEN DT5742 plugin (CAENDT5742RawEvent2StdEventConverter.cc)
        // encodes each pixel's raw hardware channel into its aux-info
        // string as "<dutname>:CH<N>:col<C>:row<R>" - there's no separate
        // structured accessor for it, so parse it out of that string.
        // Returns -1 if the pixel has no aux info, or it doesn't contain a
        // "CH<digits>" token (every non-CAEN sensor).
        auto extract_channel = [](const eudaq::StandardPlane& plane, size_t p) -> int32_t {
            if (!plane.HasPixelAuxInfo(static_cast<uint32_t>(p))) return -1;
            std::string aux = plane.GetPixelAuxInfo(static_cast<uint32_t>(p));
            auto pos = aux.find("CH");
            if (pos == std::string::npos) return -1;
            pos += 2;
            size_t end = pos;
            while (end < aux.size() && std::isdigit(static_cast<unsigned char>(aux[end]))) ++end;
            if (end == pos) return -1;
            return std::stoi(aux.substr(pos, end - pos));
        };

        // fallback_ts is the window's resolved authoritative (TLU-clock)
        // timestamp (see close_window() below) - NOT recomputed per
        // sub-event here, because a given sub-event's own GetTimeBegin()/
        // GetTimeEnd() is only ever nonzero for TLU, which contributes no
        // planes/hits of its own.
        auto append_hits = [&](const AtomicEvent& a, double fallback_ts) {
            for (size_t i = 0; i < a.evt->NumPlanes(); ++i) {
                const auto& plane = a.evt->GetPlane(i);
                std::string s_name = plane.Sensor();
                s_name.erase(std::remove_if(s_name.begin(), s_name.end(), ::isspace), s_name.end());
                if (sensor_dict_.find(s_name) == sensor_dict_.end()) sensor_dict_[s_name] = next_sensor_id_++;

                uint32_t current_sid = sensor_dict_[s_name];
                uint32_t current_pid = plane.ID();

                // MIMOSA26 (NI producer) and RD53B/croc (CMSIT producer)
                // never provide their own per-pixel timestamp
                // (plane.GetTimestamp() is always 0 for both). corry's own
                // EventLoaderEUDAQ2 runs against both with
                // `sync_by_event = true` in this project's reference
                // configs, which per is_within_event()/
                // get_trigger_position() makes corry discard TLU's hardware
                // clock for these two producers and instead assign each
                // pixel a synthetic per-trigger timestamp derived from the
                // event number: (event_number - 0.2) / 1000 [ns] +
                // detector.timeOffset() (0 here). Both producers must use
                // this same synthetic clock to match corry's output and to
                // keep RD53B/MIMOSA timestamp comparisons meaningful.
                bool uses_sync_by_event_timestamp = (s_name == "MIMOSA26" || s_name == "RD53B");

                for (size_t p = 0; p < plane.HitPixels(); ++p) {
                    // Raw EUDAQ hardware trigger number (same value already
                    // stored in /events/event_number), not the sequential
                    // global_id below - matches corry's own event numbering
                    // (which is this same hardware counter, gaps and all),
                    // so cross-referencing a specific event between corry's
                    // output and ours doesn't require guessing at an offset.
                    // global_id itself stays a purely internal loop counter
                    // (flush cadence, max_events, progress).
                    b_hi_evid.push_back(a.event_number);
                    b_hi_sens.push_back(current_sid);
                    b_hi_plane.push_back(current_pid);
                    b_hi_col.push_back(plane.GetX(p));
                    b_hi_row.push_back(plane.GetY(p));
                    int32_t channel = extract_channel(plane, p);
                    b_hi_channel.push_back(channel);
                    double ts;
                    if (plane.GetTimestamp(p) != 0) {
                        ts = static_cast<double>(plane.GetTimestamp(p)) / 1000.0;
                    } else if (uses_sync_by_event_timestamp) {
                        ts = (static_cast<double>(a.event_number) - 0.2) / 1000.0;
                    } else {
                        ts = fallback_ts;
                    }
                    b_hi_ts.push_back(ts);

                    bool has_wf = plane.HasWaveform(p);
                    b_hi_haswf.push_back(has_wf ? 1 : 0);
                    diag_hits++;
                    if (has_wf) diag_wf_hits++;
                    if (has_wf) {
                        std::vector<double> wf = plane.GetWaveform(p);
                        b_wf_t0.push_back(plane.GetWaveformX0(p));
                        b_wf_dt.push_back(plane.GetWaveformDX(p));
                        // Invert eudaq's own ADC12ToVolts (see kCaenAdcMid
                        // comment above) to recover the DRS4-corrected
                        // ADC-code value GetWaveform() calibrated FROM,
                        // rounded to the nearest integer so ds_wf_samples'
                        // double->int16 write conversion is an exact
                        // truncation, not a rounding-mode-dependent one.
                        double dc_offset_volts = caenDcOffsetVolts(channel);
                        double lsb = kCaenVpp / kCaenAdcLsbDen;
                        // No explicit reserve() here: calling
                        // reserve(size()+1024) before every waveform hit
                        // grows capacity with no slack, forcing a full
                        // reallocation on every single hit (O(n^2) over the
                        // buffer's lifetime). Plain push_back's amortized
                        // geometric growth needs only O(log n) reallocations
                        // total.
                        for (size_t s = 0; s < 1024; ++s) {
                            double v = s < wf.size() ? wf[s] : 0.0;
                            double adc_code = (v - dc_offset_volts) / lsb + kCaenAdcMid;
                            b_wf_samples.push_back(std::round(adc_code));
                        }
                    }
                }
            }
        };

        auto flush_buffers = [&]() {
            ds_ev_id.append(b_ev_id); b_ev_id.clear();
            ds_ev_num.append(b_ev_num); b_ev_num.clear();
            ds_hi_evid.append(b_hi_evid); b_hi_evid.clear();
            ds_hi_col.append(b_hi_col); b_hi_col.clear();
            ds_hi_row.append(b_hi_row); b_hi_row.clear();
            ds_hi_channel.append(b_hi_channel); b_hi_channel.clear();
            ds_hi_sens.append(b_hi_sens); b_hi_sens.clear();
            ds_hi_plane.append(b_hi_plane); b_hi_plane.clear();
            ds_hi_ts.append(b_hi_ts); b_hi_ts.clear();
            ds_hi_haswf.append(b_hi_haswf); b_hi_haswf.clear();
            ds_wf_t0.append(b_wf_t0); b_wf_t0.clear();
            ds_wf_dt.append(b_wf_dt); b_wf_dt.clear();
            ds_wf_samples.append2D(b_wf_samples, 1024); b_wf_samples.clear();
            h5file.flush(H5F_SCOPE_GLOBAL);
        };

        // With cfg_.master_description empty: windows form by simple
        // symmetric equality - whichever (sub-)event's number is seen next
        // both opens and closes windows, with no producer treated
        // specially. This matches corry's own EventLoaderEUDAQ2::
        // is_within_event() with sync_by_event=true only when every
        // producer's number sequence already agrees exactly. It does not
        // match corry in general: corry does not treat producers
        // symmetrically, since whichever detector loader is declared first
        // in the .conf alone defines each corryvreckan Event's boundary
        // from its own event number, and every other producer's
        // (sub-)event is classified before/during/after that boundary and
        // kept only if during. The purely symmetric approach diverges from
        // corry's grouping with a growing offset over the file, consistent
        // with MIMOSA26's rolling-shutter readout not being lockstep with
        // the fast digitizers' trigger-period cadence.
        //
        // With cfg_.master_description set, windows are instead anchored
        // purely on (sub-)events whose raw description matches it (see
        // AtomicEvent::is_master, set in decode_one above); every other
        // producer's (sub-)event is matched against whatever window is
        // currently open by exact event-number equality (what corry's own
        // 0.4-wide tolerance window reduces to for integer event numbers),
        // buffered in `lookahead` if it arrives before a matching window
        // has opened yet (the master producer can lag behind faster ones),
        // and dropped once its number falls behind the open window (a
        // stale straggler corry would never revisit either).
        bool have_window = false;
        uint32_t window_event_number = 0;
        bool using_master_anchor = !cfg_.master_description.empty();

        // A window's sub-events are buffered here instead of being appended
        // immediately as they arrive, because the one real timestamp for the
        // whole physical trigger comes from TLU's contribution specifically
        // (see decode_one above) - and TLU's position within a window's
        // sub-events isn't something this code should assume (first, last,
        // or anywhere else depends on how the DAQ merged producers). Buffering
        // and resolving the timestamp once the window is known to be
        // complete makes the result independent of that arrival order.
        std::vector<AtomicEvent> window_events;

        // Non-master (sub-)events seen ahead of the currently open window,
        // waiting for a matching master (sub-)event to open their window -
        // only used when using_master_anchor is true. Stays small in
        // practice (the per-producer number gap doesn't grow unboundedly),
        // so a plain vector rescanned on every new window is fine.
        std::vector<AtomicEvent> lookahead;

        auto close_window = [&]() {
            double window_ts = 0.0;
            for (const auto& we : window_events) {
                uint64_t tb = we.evt->GetTimeBegin(), te = we.evt->GetTimeEnd();
                if (tb != 0 || te != 0) {
                    window_ts = static_cast<double>(tb + te) / 2000.0;
                    break;
                }
            }
            for (const auto& we : window_events) append_hits(we, window_ts);
            window_events.clear();
        };

        // Bumps the window counter, prints progress/diagnostics, and
        // flushes periodically. Shared by both the master-anchored and
        // legacy symmetric windowing paths below, which call it from
        // different points. Returns true once max_events has been reached.
        auto advance_after_close = [&]() -> bool {
            global_id++;
            if (global_id % 100 == 0) {
                auto now = std::chrono::steady_clock::now();
                double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(now - start_time_).count();

                double elapsed_ms = elapsed * 1000.0;
                uint64_t d_ev = global_id - diag_prev_global_id;
                std::ostringstream diag_line;
                if (d_ev > 0) {
                    uint64_t d_hits = diag_hits - diag_prev_hits;
                    uint64_t d_wf = diag_wf_hits - diag_prev_wf;
                    uint64_t d_sub = diag_subevents_decoded - diag_prev_subevents;
                    double d_decode_ms = diag_decode_ms - diag_prev_decode_ms;
                    double d_wall_ms = elapsed_ms - diag_prev_elapsed_ms;
                    diag_line << "[diag] events " << diag_prev_global_id << "-" << global_id
                              << ": " << std::fixed << std::setprecision(2)
                              << (double)d_hits / d_ev << " hits/ev, "
                              << (double)d_wf / d_ev << " wf-hits/ev, "
                              << d_sub << " subevents decoded, "
                              << "wall " << d_wall_ms << " ms, decode " << d_decode_ms << " ms, "
                              << "other " << (d_wall_ms - d_decode_ms) << " ms";
                }
                printProgress(global_id, cfg_.max_events, elapsed, diag_line.str());

                diag_prev_global_id = global_id;
                diag_prev_hits = diag_hits;
                diag_prev_wf = diag_wf_hits;
                diag_prev_subevents = diag_subevents_decoded;
                diag_prev_decode_ms = diag_decode_ms;
                diag_prev_elapsed_ms = elapsed_ms;
            }
            // Flushing calls h5file.flush(H5F_SCOPE_GLOBAL), an OS-level
            // sync of the whole file, and every flush also risks
            // re-compressing (see AppendableDataset.cpp) any dataset chunk
            // that isn't exactly full yet - both costs scale with flush
            // count, not data volume. A busy CAEN board's waveform buffer
            // (b_wf_samples) can add on the order of 100+ KiB/event, so the
            // interval can't be pushed arbitrarily high either. 10000 keeps
            // the in-RAM buffer in the hundreds-of-MB range while keeping
            // flush count low.
            if (global_id % 10000 == 0) flush_buffers();
            return cfg_.max_events > 0 && global_id >= cfg_.max_events;
        };

        auto open_window = [&](uint32_t event_number, AtomicEvent&& first) {
            window_event_number = event_number;
            have_window = true;
            // event_id is the sequential window index (gap-free, matches
            // how many windows have been closed so far); event_number is
            // the raw hardware trigger number of this window (can have
            // gaps - dropped triggers are normal at the DAQ level). Kept
            // distinct so comparing the two columns reveals dropped
            // triggers directly. Only /hits/event_id (see append_hits)
            // uses the hardware number, since that's the one
            // EventLoaderHDF5 reads into Warlock's internal event numbering.
            b_ev_id.push_back(global_id);
            b_ev_num.push_back(window_event_number);
            window_events.push_back(std::move(first));
        };

        if (using_master_anchor) {
            std::cout << "[minion] Anchoring window formation on raw producer: \""
                      << cfg_.master_description << "\"" << std::endl;
        }

        bool stop = false;
        while (!stop) {
            auto opt = next_atomic();
            if (!opt) break;
            AtomicEvent cur = std::move(*opt);

            if (!using_master_anchor) {
                // Legacy symmetric windowing: whichever producer's number
                // is seen next both opens and closes windows equally.
                if (!have_window) {
                    open_window(cur.event_number, std::move(cur));
                    continue;
                }
                if (cur.event_number == window_event_number) {
                    window_events.push_back(std::move(cur));
                    continue;
                }
                if (cur.event_number < window_event_number) {
                    continue; // stale straggler, dropped
                }
                close_window();
                stop = advance_after_close();
                have_window = false;
                if (!stop) open_window(cur.event_number, std::move(cur));
                continue;
            }

            // Master-anchored windowing.
            if (cur.is_master) {
                if (have_window && cur.event_number == window_event_number) {
                    // Same window, another (sub-)event from the master
                    // producer itself (shouldn't normally happen, but no
                    // reason to treat it differently from any other
                    // contribution if it does).
                    window_events.push_back(std::move(cur));
                    continue;
                }
                if (have_window && cur.event_number < window_event_number) {
                    // Master went backwards - shouldn't happen for a
                    // well-formed file; drop as an orphaned straggler
                    // rather than reopening an already-closed window.
                    continue;
                }
                // New master (sub-)event: close whatever was open first.
                if (have_window) {
                    close_window();
                    stop = advance_after_close();
                    have_window = false;
                    if (stop) break;
                }
                // Pull out any lookahead entries that match this new
                // window right away, and drop any now provably orphaned
                // (behind the new window - the master skipped past them
                // entirely, so no window will ever claim them). Anything
                // still ahead of the new window stays buffered.
                std::vector<AtomicEvent> matched, still_pending;
                still_pending.reserve(lookahead.size());
                for (auto& le : lookahead) {
                    if (le.event_number == cur.event_number) matched.push_back(std::move(le));
                    else if (le.event_number > cur.event_number) still_pending.push_back(std::move(le));
                    // else: orphaned straggler, dropped
                }
                lookahead = std::move(still_pending);
                open_window(cur.event_number, std::move(cur));
                for (auto& m : matched) window_events.push_back(std::move(m));
                continue;
            }

            // Non-master contribution.
            if (have_window && cur.event_number == window_event_number) {
                window_events.push_back(std::move(cur));
            } else if (have_window && cur.event_number < window_event_number) {
                // Stale straggler for an already-closed window - dropped,
                // matching corry's BEFORE classification.
            } else {
                // Ahead of the currently open window (or no window open
                // yet) - the master producer hasn't caught up. Buffer
                // until it does.
                lookahead.push_back(std::move(cur));
            }
        }

        // The last window in the file never sees a following higher event
        // number to trigger close_window() above, so its buffered sub-events
        // are still sitting in window_events - close it explicitly here
        // before the final flush.
        close_window();
        flush_buffers();

        H5::Exception::dontPrint();
        for (const auto& [name, id] : sensor_dict_) {
            std::string attr_name = "dict_" + std::to_string(id);
            if (grp_hi.attrExists(attr_name)) continue;
            H5::StrType st(H5::PredType::C_S1, name.length() + 1);
            H5::Attribute attr = grp_hi.createAttribute(attr_name, st, H5S_SCALAR);
            attr.write(st, name);
        }
        // One per CAEN channel that actually produced a waveform this run -
        // the resolved (tag-or-default) Volts value ADC->Volts calibration
        // was inverted with when `samples` was encoded (see caenDcOffsetVolts
        // above), persisted so WaveformProcessingCAEN can reverse the exact same
        // calibration without re-parsing raw-file tags itself.
        for (const auto& [channel, dc_offset_volts] : caen_channel_dc_offset_volts) {
            std::string attr_name = "caen_dc_offset_volts_ch" + std::to_string(channel);
            if (grp_hi.attrExists(attr_name)) continue;
            H5::Attribute attr = grp_hi.createAttribute(attr_name, H5::PredType::NATIVE_DOUBLE, H5S_SCALAR);
            attr.write(H5::PredType::NATIVE_DOUBLE, &dc_offset_volts);
        }

        std::cout << "\n[minion] Finished! Total triggers: " << (global_id + (have_window ? 1 : 0)) << std::endl;
    }
}
