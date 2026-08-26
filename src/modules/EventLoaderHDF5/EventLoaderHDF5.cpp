#include "EventLoaderHDF5.hpp"
#include "core/ModuleFactory.hpp"
#include "core/utils/Logger.hpp"
#include "core/detector/GeometryManager.hpp"
#include <filesystem>
#include <algorithm>
#include <Eigen/Dense>
#include <chrono>
#include <sstream>
#include <future>
#include <map>

namespace framework {

    namespace {
        // Cumulative timing split between the HDF5 hyperslab read
        // (I/O-bound) and the per-hit Pixel/Waveform construction loop
        // (CPU-bound; parallelized in pass 2 of run() - see there for why
        // the compact-format row cursor must stay sequential in pass 1).
        double stat_time_read_s = 0.0;
        double stat_time_build_s = 0.0;

        // Stable sentinel a HitPlan::name can point to when no route
        // matched this hit's raw identity/row - never written to, so safe
        // to share across every hit and every batch.
        const std::string kUnknownDetectorName;
    }

    EventLoaderHDF5::EventLoaderHDF5(Configuration cfg, Configuration g_cfg, ThreadPool* pool) 
        : Module(std::move(cfg), std::move(g_cfg), pool) {
        
        file_names_ = config.getArray<std::string>("file_name");
        if (file_names_.empty()) file_names_.push_back(config.get<std::string>("file_name", "output.h5"));
        
        max_events_ = global_config.getArray<size_t>("number_of_events");
        if (max_events_.empty()) max_events_.push_back(global_config.get<size_t>("number_of_events", 0));
        
        while (max_events_.size() < file_names_.size()) max_events_.push_back(max_events_.back());

        load_waveforms_ = config.get<bool>("load_waveforms", false);
        type_filter_ = config.get<std::string>("type", "");
        // Accepts either detector names or detector types, as a single
        // bare value or a list - same convention as Tracking4D::
        // exclude_detectors_/AlignmentMillepede::exclude_detectors_.
        // Independent of type_filter_ (both apply if both are set): lets a
        // stage that needs several types loaded together (e.g. "mimosa26"
        // + "cmsit", which type_filter_ alone can't express as a single
        // string) still leave one type out, e.g.
        // exclude_detectors = "caendt5742".
        exclude_detectors_ = config.getArray<std::string>("exclude_detectors");
    }

    bool EventLoaderHDF5::isExcluded(const DetectorGeo& det) const {
        return std::find(exclude_detectors_.begin(), exclude_detectors_.end(), det.name) != exclude_detectors_.end() ||
               std::find(exclude_detectors_.begin(), exclude_detectors_.end(), det.type) != exclude_detectors_.end();
    }

    void EventLoaderHDF5::initialize() {
        if (file_names_.empty()) throw std::runtime_error("No input files specified in configuration.");
        buildRouteTable();
        openNextFile();
    }

    void EventLoaderHDF5::buildRouteTable() {
        auto& geo = GeometryManager::getInstance();
        for (auto const& name : geo.getDetectorNames()) {
            auto& det = geo.getDetector(name);
            const std::string& raw_identity = det.source_sensor.empty() ? name : det.source_sensor;
            if ((!type_filter_.empty() && det.type != type_filter_) || isExcluded(det)) {
                // Present in the .geo file, just excluded by this loader's
                // own `type`/`exclude_detectors` - track it so the per-hit
                // warning below can tell that apart from a genuinely
                // unknown sensor identity (one absent from the .geo file
                // entirely), which is the only case that actually
                // indicates a config problem.
                filtered_by_type_.insert(raw_identity);
                continue;
            }
            if (!det.source_sensor.empty()) {
                route_table_[det.source_sensor].push_back({name, det.source_channels});
            } else {
                route_table_[name].push_back({name, {}});
            }
        }
    }

    bool EventLoaderHDF5::openNextFile() {
        if (current_file_idx_ >= file_names_.size()) return false;
        
        std::string fn = file_names_[current_file_idx_];
        if (!std::filesystem::exists(fn)) throw std::runtime_error("HDF5 File not found: " + fn);
        
        try {
            H5::Exception::dontPrint();
            current_file_ = H5::H5File(fn, H5F_ACC_RDONLY);
            H5::Group hits_grp = current_file_.openGroup("/hits");

            sensor_map_.clear();
            channel_dc_offset_volts_.clear();
            static const std::string kDcOffsetPrefix = "caen_dc_offset_volts_ch";
            for (int i = 0; i < hits_grp.getNumAttrs(); ++i) {
                H5::Attribute attr = hits_grp.openAttribute(i);
                const std::string attr_name = attr.getName();
                if (attr_name.find("dict_") == 0) {
                    uint32_t id = std::stoul(attr_name.substr(5));
                    std::string sensor_type;
                    attr.read(attr.getStrType(), sensor_type);
                    sensor_map_[id] = sensor_type;
                } else if (attr_name.find(kDcOffsetPrefix) == 0) {
                    int channel = std::stoi(attr_name.substr(kDcOffsetPrefix.size()));
                    double dc_offset_volts = 0.0;
                    attr.read(H5::PredType::NATIVE_DOUBLE, &dc_offset_volts);
                    channel_dc_offset_volts_[channel] = dc_offset_volts;
                }
            }

            H5::DataSet ds_ev_id = current_file_.openDataSet("/hits/event_id");
            hsize_t dims[1];
            ds_ev_id.getSpace().getSimpleExtentDims(dims);
            total_hits_ = dims[0];
            current_hit_idx_ = 0;

            auto safe_read = [&](const std::string& path, auto& vec, const H5::DataType& type) {
                try {
                    H5::DataSet ds = current_file_.openDataSet(path);
                    vec.resize(total_hits_);
                    ds.read(vec.data(), type);
                } catch (...) {
                    throw std::runtime_error("CRITICAL: Dataset '" + path + "' missing from HDF5!");
                }
            };

            safe_read("/hits/event_id",  h_event_ids_,  H5::PredType::NATIVE_UINT64);
            // Shift this file's raw event_id numbering (which independently
            // starts near 0, same as every other file) into its own
            // reserved slice of the global range - see
            // current_file_event_offset_'s docs for why this matters once
            // more than one file is loaded.
            current_file_event_offset_ = next_file_event_offset_;
            if (current_file_event_offset_ > 0) {
                for (auto& id : h_event_ids_) id += current_file_event_offset_;
            }
            // h_event_ids_ is already globally offset at this point, so its
            // own max (not this file's local max) is the next file's base.
            next_file_event_offset_ = h_event_ids_.empty() ? current_file_event_offset_
                : *std::max_element(h_event_ids_.begin(), h_event_ids_.end()) + 1;
            safe_read("/hits/col",       h_cols_,       H5::PredType::NATIVE_INT32);
            safe_read("/hits/row",       h_rows_,       H5::PredType::NATIVE_INT32);
            safe_read("/hits/sensor_id", h_sensor_ids_, H5::PredType::NATIVE_UINT32);
            safe_read("/hits/plane_id",  h_plane_ids_,  H5::PredType::NATIVE_UINT32);
            safe_read("/hits/timestamp", h_timestamps_, H5::PredType::NATIVE_DOUBLE);

            // Optional: only present in files converted after channel-based
            // detector splitting was added. Older files fall back to -1 for
            // every hit, matching what minionRawConverter itself writes for
            // any non-CAEN sensor - a channel-filtered DetectorRoute simply
            // never matches on such a file, same as if it had no channel
            // data at all.
            try {
                H5::DataSet ds = current_file_.openDataSet("/hits/channel");
                h_channels_.resize(total_hits_);
                ds.read(h_channels_.data(), H5::PredType::NATIVE_INT32);
            } catch (...) {
                h_channels_.assign(total_hits_, -1);
            }

            compact_waveform_format_ = false;
            h_has_waveform_.clear();
            h_wf_t0_compact_.clear();
            h_wf_dt_compact_.clear();
            waveform_rows_read_ = 0;

            if (load_waveforms_) {
                try {
                    ds_samples_ = current_file_.openDataSet("/hits/samples");

                    // Newer converter output stores waveform data compactly
                    // (one row per hit that has one, not one per hit
                    // overall) and marks which hits those are via this flag.
                    // Older files won't have it - fall back to the
                    // full-length, 1:1 index-aligned layout in that case.
                    try {
                        H5::DataSet ds_haswf = current_file_.openDataSet("/hits/has_waveform");
                        h_has_waveform_.resize(total_hits_);
                        ds_haswf.read(h_has_waveform_.data(), H5::PredType::NATIVE_UINT8);
                        compact_waveform_format_ = true;

                        size_t n_wf = 0;
                        for (auto v : h_has_waveform_) if (v) ++n_wf;
                        h_wf_t0_compact_.resize(n_wf);
                        h_wf_dt_compact_.resize(n_wf);
                        if (n_wf > 0) {
                            current_file_.openDataSet("/hits/wf_t0").read(h_wf_t0_compact_.data(), H5::PredType::NATIVE_DOUBLE);
                            current_file_.openDataSet("/hits/wf_dt").read(h_wf_dt_compact_.data(), H5::PredType::NATIVE_DOUBLE);
                        }
                    } catch (...) {
                        compact_waveform_format_ = false;
                        h_has_waveform_.clear();
                        safe_read("/hits/wf_t0", h_wf_t0_, H5::PredType::NATIVE_DOUBLE);
                        safe_read("/hits/wf_dt", h_wf_dt_, H5::PredType::NATIVE_DOUBLE);
                    }

                    has_waveforms_ = true;
                    WR_LOG(STATUS, std::string("Waveform analysis enabled. 2D Hyperslab reading active (") +
                                       (compact_waveform_format_ ? "compact format)." : "legacy 1:1 format)."));
                } catch (...) {
                    has_waveforms_ = false;
                    WR_LOG(WARNING, "Waveforms requested but '/hits/samples' not found in file.");
                }
            } else {
                has_waveforms_ = false;
                WR_LOG(STATUS, "Waveform loading disabled. Disk I/O optimized.");
            }

            WR_LOG(STATUS, "Loaded file " + std::to_string(current_file_idx_+1) + "/" + std::to_string(file_names_.size()) + ": " + fn + " (" + std::to_string(total_hits_) + " hits)");
            return true;

        } catch (const H5::Exception& e) {
            throw std::runtime_error("HDF5 Error in " + getName() + ": " + e.getDetailMsg());
        }
    }

    void EventLoaderHDF5::run(DataBatch& batch) {
        size_t b_size = global_config.get<size_t>("batch_size", 5000);
        // FrameworkManager shrinks this as a `number_of_tracks` target gets
        // close, so the final batch doesn't read a full extra chunk of
        // events beyond what's actually needed to reach it.
        if (batch.requested_batch_size > 0) b_size = std::min(b_size, batch.requested_batch_size);

        // Crossing a file boundary (either this file's own data is
        // exhausted, or its per-file `number_of_events` cap was just hit)
        // must NOT return an empty batch to FrameworkManager - an
        // all-empty batch is its own signal for "truly no more data
        // anywhere" (see run()'s end_of_data check), so returning early
        // here would make it stop after file 1 without ever reading
        // files 2..N, even though openNextFile() already loaded them.
        // Loop internally until there's real data to build a batch from,
        // or every file is genuinely exhausted.
        uint64_t start_event_id = 0, end_event_id = 0;
        size_t batch_start_idx = 0, batch_end_idx = 0, hits_in_batch = 0;
        while (true) {
            if (current_file_idx_ >= file_names_.size()) return;

            if (current_hit_idx_ >= total_hits_) {
                current_file_idx_++;
                if (!openNextFile()) return;
                continue;
            }

            // Re-read every call (not hoisted above the loop): must match
            // whichever file current_file_idx_ points to *after* any
            // advance above, not the file this call started with.
            size_t max_ev = max_events_[current_file_idx_];

            start_event_id = h_event_ids_[current_hit_idx_];
            end_event_id = start_event_id + b_size;

            // number_of_events is configured per-file ("stop after N events of
            // THIS file"), but start/end_event_id are now in the global,
            // cross-file-offset numbering (see current_file_event_offset_) -
            // translate the local cap into this file's own slice of that
            // range before comparing.
            if (max_ev > 0) {
                uint64_t global_max_ev = current_file_event_offset_ + max_ev;
                if (end_event_id > global_max_ev) end_event_id = global_max_ev;
            }

            batch_start_idx = current_hit_idx_;
            batch_end_idx = current_hit_idx_;
            while (batch_end_idx < total_hits_ && h_event_ids_[batch_end_idx] < end_event_id) {
                batch_end_idx++;
            }

            hits_in_batch = batch_end_idx - batch_start_idx;

            if (hits_in_batch == 0) {
                current_file_idx_++;
                if (!openNextFile()) return;
                continue;
            }

            break;
        }

        // Actual distinct events covered by this call - may be less than
        // the configured/requested batch size near a file boundary, an
        // event cap, or a shrunk final batch (see requested_batch_size).
        batch.events_in_batch = end_event_id - start_event_id;

        auto t_read_start = std::chrono::steady_clock::now();
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> sample_block;
        size_t wf_batch_start = waveform_rows_read_;
        size_t wf_count_in_batch = 0;
        if (has_waveforms_) {
            if (compact_waveform_format_) {
                // Compact layout: only hits with has_waveform==true have a
                // row at all, in the same relative order they appear among
                // hits - count how many of those fall in this batch and
                // pull exactly that many rows starting from the running
                // offset, instead of one row per hit regardless of type.
                for (size_t i = 0; i < hits_in_batch; ++i) {
                    if (h_has_waveform_[batch_start_idx + i]) ++wf_count_in_batch;
                }
                if (wf_count_in_batch > 0) {
                    sample_block.resize(wf_count_in_batch, 1024);
                    hsize_t offset[2] = { static_cast<hsize_t>(wf_batch_start), 0 };
                    hsize_t count[2]  = { static_cast<hsize_t>(wf_count_in_batch), 1024 };
                    H5::DataSpace filespace = ds_samples_.getSpace();
                    filespace.selectHyperslab(H5S_SELECT_SET, count, offset);
                    H5::DataSpace memspace(2, count);
                    ds_samples_.read(sample_block.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace);
                }
                waveform_rows_read_ += wf_count_in_batch;
            } else {
                sample_block.resize(hits_in_batch, 1024);
                hsize_t offset[2] = { static_cast<hsize_t>(batch_start_idx), 0 };
                hsize_t count[2]  = { static_cast<hsize_t>(hits_in_batch), 1024 };
                H5::DataSpace filespace = ds_samples_.getSpace();
                filespace.selectHyperslab(H5S_SELECT_SET, count, offset);
                H5::DataSpace memspace(2, count);
                ds_samples_.read(sample_block.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace);
            }
        }

        stat_time_read_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_read_start).count();

        auto t_build_start = std::chrono::steady_clock::now();

        // Pass 1 (serial, cheap): resolve each hit's detector name and
        // compact-format waveform row index. Must stay single-threaded -
        // wf_row_cursor has to advance strictly in hit order, since
        // minionRawConverter writes /hits/samples indexed by a hit's
        // position among ALL waveform-carrying hits in the batch (any
        // detector, known or not). route_cache_/sensor_map_/unknown_names_
        // are also only ever written here, never in pass 2 below, so pass 2
        // can read them without synchronization.
        struct HitPlan {
            const std::string* name;
            int row_offset; // subtracted from this hit's raw row (0 unless split)
            bool has_wf;
            size_t row_idx;
        };
        std::vector<HitPlan> plan(hits_in_batch);
        size_t wf_row_cursor = 0; // position within this batch's compact sample_block
        for (size_t i = 0; i < hits_in_batch; ++i) {
            size_t idx = batch_start_idx + i;
            uint64_t cache_key = (static_cast<uint64_t>(h_sensor_ids_[idx]) << 32) | h_plane_ids_[idx];

            auto cached = route_cache_.find(cache_key);
            if (cached == route_cache_.end()) {
                std::string candidate = sensor_map_[h_sensor_ids_[idx]] + "_" + std::to_string(h_plane_ids_[idx]);
                auto route_it = route_table_.find(candidate);
                const std::vector<DetectorRoute>* routes = (route_it != route_table_.end()) ? &route_it->second : nullptr;
                cached = route_cache_.emplace(cache_key, routes).first;
                if (!routes && unknown_names_.find(candidate) == unknown_names_.end()) {
                    if (filtered_by_type_.count(candidate)) {
                        WR_LOG(STATUS, "Sensor '" + candidate + "' excluded by this loader's `type`/`exclude_detectors` config - discarding its hits.");
                    } else {
                        WR_LOG(WARNING, "HDF5 data contains sensor '" + candidate + "', but it is NOT in your .geo file! Discarding its hits.");
                    }
                    unknown_names_.insert(candidate);
                }
            }

            // A raw identity normally routes to exactly one detector (an
            // empty filter_channels, matching every hit); a split raw
            // sensor routes to whichever of its detectors' filter_channels
            // contains this specific hit's raw channel. Every channel in a
            // route's filter list is assumed to share one raw row, so
            // matching by channel rather than row also correctly separates
            // a row that itself mixes DAQ trigger groups - a plain
            // row-based filter couldn't do that.
            const std::string* resolved_name = &kUnknownDetectorName;
            int row_offset = 0;
            if (cached->second) {
                for (auto const& route : *cached->second) {
                    if (route.filter_channels.empty()) {
                        resolved_name = &route.name;
                        row_offset = 0;
                        break;
                    }
                    if (std::find(route.filter_channels.begin(), route.filter_channels.end(), h_channels_[idx]) !=
                        route.filter_channels.end()) {
                        resolved_name = &route.name;
                        row_offset = h_rows_[idx]; // collapse this die's raw row to local row 0
                        break;
                    }
                }
            }

            const std::string& type = sensor_map_[h_sensor_ids_[idx]];
            bool this_hit_has_wf = false;
            if (has_waveforms_) {
                if (compact_waveform_format_) {
                    this_hit_has_wf = h_has_waveform_[idx] != 0;
                } else if (type == "caendt5742") {
                    this_hit_has_wf = true;
                } else if (warned_unknown_wf_types_.insert(type).second) {
                    // Legacy (pre-compact) files carry no per-hit
                    // has_waveform flag, so waveform association here
                    // falls back to a hardcoded sensor-type check - an
                    // unrecognized type just means no waveform is
                    // attached for it (never garbage/misaligned data),
                    // but that's silent otherwise. Warn once per type
                    // seen, not once per hit (this runs in a per-hit hot
                    // loop).
                    WR_LOG(WARNING, "EventLoaderHDF5: sensor type '" + type + "' isn't recognized as a "
                                     "known waveform-carrying digitizer in this file's legacy "
                                     "(pre-compact) waveform format - no waveforms will be attached for "
                                     "it. If this is a real waveform digitizer, add its sensor-type "
                                     "string here.");
                }
            }

            size_t row_idx = compact_waveform_format_ ? wf_row_cursor : i;
            if (this_hit_has_wf && compact_waveform_format_) ++wf_row_cursor;

            plan[i] = { resolved_name, row_offset, this_hit_has_wf, row_idx };
        }

        // Pass 2 (parallel): the actual heap allocation / waveform-sample
        // copy work, now safe to spread across threads since every hit's
        // name/row_idx was already resolved sequentially above and nothing
        // here mutates shared state outside its own thread-local buffers
        // until the final merge under batch_mtx_.
        size_t n_threads = thread_pool->getThreadCount();
        size_t chunk_size = (hits_in_batch + n_threads - 1) / n_threads;
        if (chunk_size == 0) chunk_size = 1;
        std::vector<std::future<void>> build_futures;

        for (size_t chunk_start = 0; chunk_start < hits_in_batch; chunk_start += chunk_size) {
            size_t chunk_end = std::min(chunk_start + chunk_size, hits_in_batch);

            build_futures.push_back(thread_pool->submit([this, chunk_start, chunk_end, batch_start_idx,
                                                           wf_batch_start, &plan, &sample_block, &batch]() {
                std::map<std::string, std::vector<std::shared_ptr<Pixel>>> local_pixels;
                std::map<std::string, std::vector<std::shared_ptr<Waveform>>> local_waveforms;

                for (size_t i = chunk_start; i < chunk_end; ++i) {
                    if (plan[i].name->empty()) continue;
                    size_t idx = batch_start_idx + i;
                    double ts = h_timestamps_[idx];
                    uint64_t ev_id = h_event_ids_[idx];
                    const std::string& name = *plan[i].name;
                    // Raw row remapped into this detector's own local
                    // numbering - identity (row_offset 0) unless `name` is
                    // one detector of a raw sensor split by channel (see
                    // DetectorGeo::source_sensor/source_channels), in which
                    // case every hit routed here already shares one raw row
                    // (every channel in a route's filter list is on the
                    // same physical die), so this always lands at local row 0.
                    int row = h_rows_[idx] - plan[i].row_offset;
                    int32_t channel = h_channels_[idx];

                    if (plan[i].has_wf) {
                        size_t row_idx = plan[i].row_idx;
                        std::vector<double> samples(1024);
                        Eigen::VectorXd wf_row = sample_block.row(row_idx);
                        std::copy(wf_row.data(), wf_row.data() + 1024, samples.begin());

                        auto wf = std::make_shared<Waveform>(name, h_cols_[idx], row, ts, std::move(samples), ev_id, channel);
                        auto dc_it = channel_dc_offset_volts_.find(channel);
                        if (dc_it != channel_dc_offset_volts_.end()) wf->dc_offset_volts = dc_it->second;
                        if (compact_waveform_format_) {
                            wf->t0 = h_wf_t0_compact_[wf_batch_start + row_idx];
                            wf->dt = h_wf_dt_compact_[wf_batch_start + row_idx];
                        } else {
                            wf->t0 = h_wf_t0_[idx];
                            wf->dt = h_wf_dt_[idx];
                        }
                        local_waveforms[name].push_back(std::move(wf));
                    } else {
                        // Matches corry's EventLoaderEUDAQ2, which drops a
                        // masked pixel's hits at load time rather than
                        // letting them reach the clipboard - every
                        // downstream module sees an already-clean pixel
                        // stream. ClusteringSpatial re-checks masking on its
                        // own accord, but any module reading batch.pixels
                        // directly (e.g. Correlations' hitmap) depends on
                        // this filter to match corry's own output, which
                        // never includes masked/noisy pixels' hits.
                        auto& det = GeometryManager::getInstance().getDetector(name);
                        if (det.isMasked(h_cols_[idx], row)) continue;
                        local_pixels[name].push_back(std::make_shared<Pixel>(name, h_cols_[idx], row, 1.0, ts, ev_id, channel));
                    }
                }

                std::lock_guard<std::mutex> lock(batch_mtx_);
                for (auto& [name, vec] : local_pixels) {
                    batch.pixels[name].insert(batch.pixels[name].end(), vec.begin(), vec.end());
                }
                for (auto& [name, vec] : local_waveforms) {
                    batch.waveforms[name].insert(batch.waveforms[name].end(), vec.begin(), vec.end());
                }
            }));
        }
        for (auto& f : build_futures) f.get();

        stat_time_build_s += std::chrono::duration<double>(std::chrono::steady_clock::now() - t_build_start).count();
        current_hit_idx_ = batch_end_idx;
    }

    void EventLoaderHDF5::finalize() {
        std::ostringstream stats;
        stats << "EventLoaderHDF5 phase timing (seconds): "
              << "hyperslab_read=" << stat_time_read_s
              << " object_build=" << stat_time_build_s;
        WR_LOG(STATUS, stats.str());
    }

    REGISTER_MODULE(EventLoaderHDF5)
}