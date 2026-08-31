/** @file Reader.cpp */
#include "Reader.hpp"
#include "core/ModuleFactory.hpp"
#include "core/utils/Logger.hpp"
#include "objects/Track.hpp"
#include "objects/Cluster.hpp"
#include "objects/Waveform.hpp"
#include "core/detector/GeometryManager.hpp"
#include <algorithm>
#include <cmath>

namespace framework {

    namespace {
        /// Reads dataset `name` from `grp` in full, or returns an empty
        /// vector if it doesn't exist - used only for the cheap per-row
        /// *index* columns (event_id/track_id/cluster_id) that are read
        /// eagerly; the heavy payload columns are never read this way (see
        /// gatherRows() below).
        template <typename T>
        std::vector<T> readVec(H5::Group& grp, const std::string& name, H5::PredType type) {
            if (!grp.nameExists(name)) return {};
            H5::DataSet ds = grp.openDataSet(name);
            H5::DataSpace space = ds.getSpace();
            hsize_t dims[1] = {0};
            space.getSimpleExtentDims(dims);
            std::vector<T> data(dims[0]);
            if (dims[0] > 0) ds.read(data.data(), type);
            return data;
        }

        /// Reads a 1D numeric attribute off `grp`, or an empty vector if
        /// `name` isn't present.
        std::vector<double> readAttrVec(H5::Group& grp, const std::string& name) {
            if (!grp.attrExists(name)) return {};
            H5::Attribute attr = grp.openAttribute(name);
            H5::DataSpace space = attr.getSpace();
            hsize_t dims[1] = {0};
            space.getSimpleExtentDims(dims);
            std::vector<double> data(dims[0]);
            if (dims[0] > 0) attr.read(H5::PredType::NATIVE_DOUBLE, data.data());
            return data;
        }

        /// Fetches exactly the requested (possibly unsorted, possibly
        /// duplicated) rows of a rank-1 dataset via one batched hyperslab
        /// read, returned as a row_idx -> value map for the caller to look
        /// up by row. Adjacent requested rows are coalesced into single
        /// contiguous blocks before the read (rows referenced by one
        /// Reader batch tend to cluster together on disk, since they were
        /// all written close together in time), so this stays efficient
        /// even though it's not a single contiguous range in general -
        /// this is the mechanism that lets Reader avoid ever materializing
        /// a whole dataset in memory, see Reader.hpp's class docs.
        template <typename T>
        std::map<uint32_t, T> gatherRows(H5::DataSet& ds, const std::vector<uint32_t>& requested_rows,
                                          H5::PredType type) {
            std::map<uint32_t, T> out;
            if (requested_rows.empty()) return out;
            std::vector<uint32_t> rows = requested_rows;
            std::sort(rows.begin(), rows.end());
            rows.erase(std::unique(rows.begin(), rows.end()), rows.end());

            H5::DataSpace filespace = ds.getSpace();
            bool first = true;
            size_t i = 0;
            while (i < rows.size()) {
                size_t j = i + 1;
                while (j < rows.size() && rows[j] == rows[j - 1] + 1) ++j;
                hsize_t offset[1] = {rows[i]};
                hsize_t count[1] = {j - i};
                filespace.selectHyperslab(first ? H5S_SELECT_SET : H5S_SELECT_OR, count, offset);
                first = false;
                i = j;
            }
            hsize_t n = rows.size();
            H5::DataSpace memspace(1, &n);
            std::vector<T> buf(n);
            ds.read(buf.data(), type, memspace, filespace);
            for (size_t k = 0; k < rows.size(); ++k) out[rows[k]] = buf[k];
            return out;
        }

        /// Same idea as gatherRows(), for a rank-2 (fixed-`width`-per-row)
        /// dataset - used for start_time_cfd.
        std::map<uint32_t, std::vector<double>> gatherRows2D(H5::DataSet& ds,
                                                               const std::vector<uint32_t>& requested_rows,
                                                               hsize_t width) {
            std::map<uint32_t, std::vector<double>> out;
            if (requested_rows.empty() || width == 0) return out;
            std::vector<uint32_t> rows = requested_rows;
            std::sort(rows.begin(), rows.end());
            rows.erase(std::unique(rows.begin(), rows.end()), rows.end());

            H5::DataSpace filespace = ds.getSpace();
            bool first = true;
            size_t i = 0;
            while (i < rows.size()) {
                size_t j = i + 1;
                while (j < rows.size() && rows[j] == rows[j - 1] + 1) ++j;
                hsize_t offset[2] = {rows[i], 0};
                hsize_t count[2] = {j - i, width};
                filespace.selectHyperslab(first ? H5S_SELECT_SET : H5S_SELECT_OR, count, offset);
                first = false;
                i = j;
            }
            hsize_t n = rows.size();
            hsize_t mem_dims[2] = {n, width};
            H5::DataSpace memspace(2, mem_dims);
            std::vector<double> buf(n * width);
            ds.read(buf.data(), H5::PredType::NATIVE_DOUBLE, memspace, filespace);
            for (size_t k = 0; k < rows.size(); ++k)
                out[rows[k]] = std::vector<double>(buf.begin() + static_cast<long>(k * width),
                                                     buf.begin() + static_cast<long>((k + 1) * width));
            return out;
        }
    }

    Reader::Reader(Configuration cfg, Configuration g_cfg, ThreadPool* pool)
        : Module(std::move(cfg), std::move(g_cfg), pool) {
        filenames_ = config.getArray<std::string>("filenames");
        if (filenames_.empty()) {
            std::string single = config.get<std::string>("filename", "");
            if (!single.empty()) filenames_.push_back(single);
        }
        batch_size_ = global_config.get<size_t>("batch_size", 5000);
    }

    void Reader::initialize() {
        WR_LOG(STATUS, "Initializing Reader...");
        track_handles_.resize(filenames_.size());
        for (size_t i = 0; i < filenames_.size(); ++i) loadFile(filenames_[i], i);

        // Index rows must come out ordered by event_id, so that batches
        // emitted one after another cover strictly increasing event ranges
        // - the same guarantee EventLoaderHDF5 provides from the raw hit
        // stream, which downstream modules (and events_in_batch
        // bookkeeping) rely on. On-disk row order can't be assumed sorted
        // (Tracking4D/WaveformProcessingCAEN don't guarantee it - see
        // AppendableDataset's writers), so this sort only touches the
        // cheap index entries, not the full row payload.
        std::sort(track_index_.begin(), track_index_.end(),
                  [](const TrackIndex& a, const TrackIndex& b) { return a.event_id < b.event_id; });

        for (auto& [det_name, rows] : waveforms_by_det_) {
            std::sort(rows.begin(), rows.end(),
                      [](const WaveformIndex& a, const WaveformIndex& b) { return a.event_id < b.event_id; });
        }

        // Builds tc_track_ids_/tc_offsets_/tc_cluster_ids_ from the raw
        // pairs accumulated across every loadFile() call above, and marks
        // ClusterLoc::linked for every cluster_id referenced by any link
        // - see buildTrackClusterIndex()'s and ClusterLoc::linked's own
        // docs for why this replaced a std::map<uint64_t,
        // std::vector<uint64_t>> plus a separate std::set<uint64_t> of
        // linked ids (both gone now).
        buildTrackClusterIndex();

        // Clusters never referenced by any /track_cluster_links entry are
        // raw DUT clusters saved by Tracking4D::saveDutClusters() (e.g.
        // RD53B_115) - group and sort them the same way as waveforms, so
        // run() can replay them event-synced, fetching full rows through
        // the same fetchClusters() path as track-linked ones.
        size_t n_clusters_total = 0;
        for (uint64_t cid = 0; cid < cluster_index_.size(); ++cid) {
            const auto& loc = cluster_index_[cid];
            if (!loc.valid) continue;
            ++n_clusters_total;
            if (loc.linked) continue;
            unlinked_clusters_by_det_[cluster_det_names_[loc.det_idx]].push_back({loc.event_id, cid});
        }
        for (auto& [det_name, rows] : unlinked_clusters_by_det_) {
            std::sort(rows.begin(), rows.end(),
                      [](const ClusterIndexEntry& a, const ClusterIndexEntry& b) { return a.event_id < b.event_id; });
        }

        // Every container above is now sorted by event_id, so its own max
        // is just its last element - overall max is the largest of those.
        // Used to give the track_index_.empty() single-shot replay path
        // (see run()) an accurate events_in_batch too.
        if (!track_index_.empty()) max_loaded_event_id_ = std::max(max_loaded_event_id_, track_index_.back().event_id);
        for (auto const& [det_name, rows] : waveforms_by_det_) {
            if (!rows.empty()) max_loaded_event_id_ = std::max(max_loaded_event_id_, rows.back().event_id);
        }
        for (auto const& [det_name, rows] : unlinked_clusters_by_det_) {
            if (!rows.empty()) max_loaded_event_id_ = std::max(max_loaded_event_id_, rows.back().event_id);
        }

        WR_LOG(STATUS, "Reader: indexed " + std::to_string(track_index_.size()) + " tracks, " +
                    std::to_string(n_clusters_total) + " clusters across " +
                    std::to_string(waveforms_by_det_.size()) + " waveform channels, " +
                    std::to_string(unlinked_clusters_by_det_.size()) +
                    " unassociated DUT cluster channels (heavy fields streamed on demand, not preloaded).");
    }

    void Reader::buildTrackClusterIndex() {
        std::sort(track_cluster_pairs_raw_.begin(), track_cluster_pairs_raw_.end());

        tc_track_ids_.reserve(track_cluster_pairs_raw_.size());
        tc_offsets_.reserve(track_cluster_pairs_raw_.size() + 1);
        tc_cluster_ids_.reserve(track_cluster_pairs_raw_.size());

        bool have_current = false;
        uint64_t current_track = 0;
        for (auto const& [track_id, cluster_id] : track_cluster_pairs_raw_) {
            if (!have_current || track_id != current_track) {
                tc_track_ids_.push_back(track_id);
                tc_offsets_.push_back(static_cast<uint32_t>(tc_cluster_ids_.size()));
                current_track = track_id;
                have_current = true;
            }
            tc_cluster_ids_.push_back(cluster_id);
            if (cluster_id < cluster_index_.size()) cluster_index_[cluster_id].linked = true;
        }
        tc_offsets_.push_back(static_cast<uint32_t>(tc_cluster_ids_.size())); // end sentinel

        tc_track_ids_.shrink_to_fit();
        tc_offsets_.shrink_to_fit();
        tc_cluster_ids_.shrink_to_fit();

        // Only needed to build the CSR form above - drop it now rather
        // than holding both representations of the same data for the
        // rest of this Reader's lifetime. swap-with-empty (not clear()),
        // since clear() alone doesn't release a vector's own capacity.
        std::vector<std::pair<uint64_t, uint64_t>>().swap(track_cluster_pairs_raw_);
    }

    std::pair<const uint64_t*, size_t> Reader::clustersForTrack(uint64_t track_id) const {
        auto it = std::lower_bound(tc_track_ids_.begin(), tc_track_ids_.end(), track_id);
        if (it == tc_track_ids_.end() || *it != track_id) return {nullptr, 0};
        size_t idx = static_cast<size_t>(it - tc_track_ids_.begin());
        return {tc_cluster_ids_.data() + tc_offsets_[idx], tc_offsets_[idx + 1] - tc_offsets_[idx]};
    }

    void Reader::loadFile(const std::string& filename, size_t file_idx) {
        files_.push_back(std::make_unique<H5::H5File>(filename, H5F_ACC_RDONLY));
        H5::H5File& file = *files_.back();

        if (file.nameExists("/tracks")) {
            H5::Group tr = file.openGroup("/tracks");
            auto event_id = readVec<uint64_t>(tr, "event_id", H5::PredType::NATIVE_UINT64);
            auto track_id = readVec<uint64_t>(tr, "track_id", H5::PredType::NATIVE_UINT64);

            track_handles_[file_idx] = TrackFileHandles{tr.openDataSet("mx"), tr.openDataSet("cx"),
                                                          tr.openDataSet("my"), tr.openDataSet("cy"),
                                                          tr.openDataSet("chi2"), tr.openDataSet("ndof")};

            for (size_t i = 0; i < event_id.size(); ++i) {
                track_index_.push_back(
                    {event_id[i], track_id[i], static_cast<uint32_t>(file_idx), static_cast<uint32_t>(i)});
            }
        }

        if (file.nameExists("/track_cluster_links")) {
            H5::Group lg = file.openGroup("/track_cluster_links");
            auto link_track = readVec<uint64_t>(lg, "track_id", H5::PredType::NATIVE_UINT64);
            auto link_cluster = readVec<uint64_t>(lg, "cluster_id", H5::PredType::NATIVE_UINT64);
            track_cluster_pairs_raw_.reserve(track_cluster_pairs_raw_.size() + link_track.size());
            for (size_t i = 0; i < link_track.size(); ++i) {
                track_cluster_pairs_raw_.emplace_back(link_track[i], link_cluster[i]);
            }
        }

        if (file.nameExists("/clusters")) {
            H5::Group clusters_grp = file.openGroup("/clusters");
            hsize_t n_det = clusters_grp.getNumObjs();
            for (hsize_t i = 0; i < n_det; ++i) {
                std::string det_name = clusters_grp.getObjnameByIdx(i);
                H5::Group dg = clusters_grp.openGroup(det_name);

                auto event_id = readVec<uint64_t>(dg, "event_id", H5::PredType::NATIVE_UINT64);
                auto cluster_id = readVec<uint64_t>(dg, "cluster_id", H5::PredType::NATIVE_UINT64);

                size_t det_idx;
                auto it = cluster_det_idx_.find(det_name);
                if (it == cluster_det_idx_.end()) {
                    det_idx = cluster_det_names_.size();
                    cluster_det_names_.push_back(det_name);
                    cluster_det_idx_[det_name] = det_idx;
                } else {
                    det_idx = it->second;
                }

                ClusterDetHandles handles;
                handles.col = dg.openDataSet("col");
                handles.row = dg.openDataSet("row");
                handles.global_x = dg.openDataSet("global_x");
                handles.global_y = dg.openDataSet("global_y");
                handles.global_z = dg.openDataSet("global_z");
                handles.charge = dg.openDataSet("charge");
                handles.timestamp = dg.openDataSet("timestamp");
                // Optional: only present in files saved after channel
                // tracking was added to the cluster schema.
                handles.has_channel = dg.nameExists("channel");
                if (handles.has_channel) handles.channel = dg.openDataSet("channel");
                cluster_handles_[{det_idx, file_idx}] = std::move(handles);

                for (size_t j = 0; j < event_id.size(); ++j) {
                    uint64_t cid = cluster_id[j];
                    if (cid >= cluster_index_.size()) cluster_index_.resize(cid + 1);
                    cluster_index_[cid] = ClusterLoc{event_id[j], static_cast<uint32_t>(file_idx),
                                                      static_cast<uint32_t>(j), static_cast<uint16_t>(det_idx), true};
                }
            }
        }

        if (file.nameExists("/waveforms")) {
            H5::Group wf_grp = file.openGroup("/waveforms");
            hsize_t n_det = wf_grp.getNumObjs();
            for (hsize_t i = 0; i < n_det; ++i) {
                std::string det_name = wf_grp.getObjnameByIdx(i);
                H5::Group dg = wf_grp.openGroup(det_name);

                auto event_id = readVec<uint64_t>(dg, "event_id", H5::PredType::NATIVE_UINT64);

                WaveformDetHandles handles;
                handles.col = dg.openDataSet("col");
                handles.row = dg.openDataSet("row");
                // Optional: only present in files saved after channel
                // tracking was added - falls back to -1 in fetchWaveforms().
                handles.has_channel = dg.nameExists("channel");
                if (handles.has_channel) handles.channel = dg.openDataSet("channel");
                handles.timestamp = dg.openDataSet("timestamp");
                handles.baseline = dg.openDataSet("baseline");
                handles.noise = dg.openDataSet("noise");
                handles.amplitude = dg.openDataSet("amplitude");
                handles.snr = dg.openDataSet("snr");
                handles.rise_time = dg.openDataSet("rise_time");
                handles.integral = dg.openDataSet("integral");
                handles.t0 = dg.openDataSet("t0");
                handles.dt = dg.openDataSet("dt");
                handles.has_tot_ton = dg.nameExists("time_over_threshold") && dg.nameExists("time_over_noise");
                if (handles.has_tot_ton) {
                    handles.time_over_threshold = dg.openDataSet("time_over_threshold");
                    handles.time_over_noise = dg.openDataSet("time_over_noise");
                }

                // Per-CFD start_time_cfd (current format) - falls back to
                // the old scalar "start_time" dataset for files saved
                // before per-CFD storage existed, treating it as a
                // single-fraction file.
                if (dg.nameExists("start_time_cfd")) {
                    handles.start_time_cfd = dg.openDataSet("start_time_cfd");
                    H5::DataSpace sp = handles.start_time_cfd.getSpace();
                    hsize_t dims[2] = {0, 0};
                    sp.getSimpleExtentDims(dims);
                    handles.cfd_width = dims[1];
                    handles.cfd_rank2 = true;
                } else if (dg.nameExists("start_time")) {
                    handles.start_time_cfd = dg.openDataSet("start_time");
                    handles.cfd_width = 1;
                    handles.cfd_rank2 = false;
                }

                if (cfd_fractions_.empty()) {
                    auto fractions = readAttrVec(dg, "cfd_fractions");
                    if (!fractions.empty()) cfd_fractions_ = std::move(fractions);
                    else if (handles.cfd_width == 1) cfd_fractions_ = {0.7}; // legacy-file default
                }

                waveform_handles_[{det_name, file_idx}] = std::move(handles);

                auto& idx = waveforms_by_det_[det_name];
                for (size_t j = 0; j < event_id.size(); ++j) {
                    idx.push_back({event_id[j], static_cast<uint32_t>(file_idx), static_cast<uint32_t>(j)});
                }
            }
        }
    }

    std::map<uint64_t, Reader::ClusterRow> Reader::fetchClusters(const std::vector<uint64_t>& cluster_ids) {
        std::map<uint64_t, ClusterRow> out;
        if (cluster_ids.empty()) return out;

        std::map<std::pair<size_t, size_t>, std::vector<uint64_t>> cids_by_group; // (det_idx, file_idx) -> cluster_ids
        for (uint64_t cid : cluster_ids) {
            if (cid >= cluster_index_.size() || !cluster_index_[cid].valid) continue;
            const auto& loc = cluster_index_[cid];
            cids_by_group[{loc.det_idx, loc.file_idx}].push_back(cid);
        }

        for (auto& [key, cids] : cids_by_group) {
            auto hit = cluster_handles_.find(key);
            if (hit == cluster_handles_.end()) continue;
            auto& h = hit->second;

            std::vector<uint32_t> rows;
            rows.reserve(cids.size());
            for (uint64_t cid : cids) rows.push_back(cluster_index_[cid].row_idx);

            auto m_col = gatherRows<double>(h.col, rows, H5::PredType::NATIVE_DOUBLE);
            auto m_row = gatherRows<double>(h.row, rows, H5::PredType::NATIVE_DOUBLE);
            auto m_gx = gatherRows<double>(h.global_x, rows, H5::PredType::NATIVE_DOUBLE);
            auto m_gy = gatherRows<double>(h.global_y, rows, H5::PredType::NATIVE_DOUBLE);
            auto m_gz = gatherRows<double>(h.global_z, rows, H5::PredType::NATIVE_DOUBLE);
            auto m_charge = gatherRows<double>(h.charge, rows, H5::PredType::NATIVE_DOUBLE);
            auto m_timestamp = gatherRows<double>(h.timestamp, rows, H5::PredType::NATIVE_DOUBLE);
            std::map<uint32_t, int32_t> m_channel;
            if (h.has_channel) m_channel = gatherRows<int32_t>(h.channel, rows, H5::PredType::NATIVE_INT32);

            for (uint64_t cid : cids) {
                uint32_t r = cluster_index_[cid].row_idx;
                out[cid] = ClusterRow{cluster_index_[cid].event_id,
                                       m_col.at(r),
                                       m_row.at(r),
                                       m_gx.at(r),
                                       m_gy.at(r),
                                       m_gz.at(r),
                                       m_charge.at(r),
                                       m_timestamp.at(r),
                                       h.has_channel ? m_channel.at(r) : -1};
            }
        }
        return out;
    }

    std::vector<Reader::WaveformRow> Reader::fetchWaveforms(const std::string& det_name,
                                                              const std::vector<WaveformIndex>& entries) {
        std::vector<WaveformRow> out(entries.size());
        if (entries.empty()) return out;

        std::map<size_t, std::vector<uint32_t>> rows_by_file;
        for (const auto& e : entries) rows_by_file[e.file_idx].push_back(e.row_idx);

        std::map<size_t, std::map<uint32_t, int32_t>> m_col, m_row, m_channel;
        std::map<size_t, std::map<uint32_t, double>> m_timestamp, m_baseline, m_noise, m_amplitude, m_snr,
            m_rise_time, m_integral, m_t0, m_dt, m_tot, m_ton;
        std::map<size_t, std::map<uint32_t, std::vector<double>>> m_start;
        std::map<size_t, bool> has_channel_by_file, has_tot_ton_by_file;

        for (auto& [file_idx, rows] : rows_by_file) {
            auto hit = waveform_handles_.find({det_name, file_idx});
            if (hit == waveform_handles_.end()) continue;
            auto& h = hit->second;

            m_col[file_idx] = gatherRows<int32_t>(h.col, rows, H5::PredType::NATIVE_INT32);
            m_row[file_idx] = gatherRows<int32_t>(h.row, rows, H5::PredType::NATIVE_INT32);
            has_channel_by_file[file_idx] = h.has_channel;
            if (h.has_channel) m_channel[file_idx] = gatherRows<int32_t>(h.channel, rows, H5::PredType::NATIVE_INT32);
            m_timestamp[file_idx] = gatherRows<double>(h.timestamp, rows, H5::PredType::NATIVE_DOUBLE);
            m_baseline[file_idx] = gatherRows<double>(h.baseline, rows, H5::PredType::NATIVE_DOUBLE);
            m_noise[file_idx] = gatherRows<double>(h.noise, rows, H5::PredType::NATIVE_DOUBLE);
            m_amplitude[file_idx] = gatherRows<double>(h.amplitude, rows, H5::PredType::NATIVE_DOUBLE);
            m_snr[file_idx] = gatherRows<double>(h.snr, rows, H5::PredType::NATIVE_DOUBLE);
            m_rise_time[file_idx] = gatherRows<double>(h.rise_time, rows, H5::PredType::NATIVE_DOUBLE);
            m_integral[file_idx] = gatherRows<double>(h.integral, rows, H5::PredType::NATIVE_DOUBLE);
            m_t0[file_idx] = gatherRows<double>(h.t0, rows, H5::PredType::NATIVE_DOUBLE);
            m_dt[file_idx] = gatherRows<double>(h.dt, rows, H5::PredType::NATIVE_DOUBLE);
            has_tot_ton_by_file[file_idx] = h.has_tot_ton;
            if (h.has_tot_ton) {
                m_tot[file_idx] = gatherRows<double>(h.time_over_threshold, rows, H5::PredType::NATIVE_DOUBLE);
                m_ton[file_idx] = gatherRows<double>(h.time_over_noise, rows, H5::PredType::NATIVE_DOUBLE);
            }
            if (h.cfd_width > 0) {
                if (h.cfd_rank2) {
                    m_start[file_idx] = gatherRows2D(h.start_time_cfd, rows, h.cfd_width);
                } else {
                    auto flat = gatherRows<double>(h.start_time_cfd, rows, H5::PredType::NATIVE_DOUBLE);
                    for (auto& [r, v] : flat) m_start[file_idx][r] = {v};
                }
            }
        }

        for (size_t k = 0; k < entries.size(); ++k) {
            const auto& e = entries[k];
            WaveformRow wr;
            wr.event_id = e.event_id;
            wr.col = m_col.at(e.file_idx).at(e.row_idx);
            wr.row = m_row.at(e.file_idx).at(e.row_idx);
            wr.channel = has_channel_by_file[e.file_idx] ? m_channel.at(e.file_idx).at(e.row_idx) : -1;
            wr.timestamp = m_timestamp.at(e.file_idx).at(e.row_idx);
            wr.baseline = m_baseline.at(e.file_idx).at(e.row_idx);
            wr.noise = m_noise.at(e.file_idx).at(e.row_idx);
            wr.amplitude = m_amplitude.at(e.file_idx).at(e.row_idx);
            wr.snr = m_snr.at(e.file_idx).at(e.row_idx);
            wr.rise_time = m_rise_time.at(e.file_idx).at(e.row_idx);
            wr.integral = m_integral.at(e.file_idx).at(e.row_idx);
            wr.t0 = m_t0.at(e.file_idx).at(e.row_idx);
            wr.dt = m_dt.at(e.file_idx).at(e.row_idx);
            if (has_tot_ton_by_file[e.file_idx]) {
                wr.time_over_threshold = m_tot.at(e.file_idx).at(e.row_idx);
                wr.time_over_noise = m_ton.at(e.file_idx).at(e.row_idx);
            }
            auto file_start_it = m_start.find(e.file_idx);
            if (file_start_it != m_start.end()) {
                auto row_start_it = file_start_it->second.find(e.row_idx);
                if (row_start_it != file_start_it->second.end()) wr.start_times = row_start_it->second;
            }
            out[k] = std::move(wr);
        }
        return out;
    }

    std::vector<Reader::TrackFields> Reader::fetchTrackFields(const std::vector<TrackIndex>& entries) {
        std::vector<TrackFields> out(entries.size());
        if (entries.empty()) return out;

        std::map<size_t, std::vector<uint32_t>> rows_by_file;
        for (const auto& e : entries) rows_by_file[e.file_idx].push_back(e.row_idx);

        std::map<size_t, std::map<uint32_t, double>> m_mx, m_cx, m_my, m_cy, m_chi2;
        std::map<size_t, std::map<uint32_t, int32_t>> m_ndof;

        for (auto& [file_idx, rows] : rows_by_file) {
            auto& h = track_handles_.at(file_idx);
            m_mx[file_idx] = gatherRows<double>(h.mx, rows, H5::PredType::NATIVE_DOUBLE);
            m_cx[file_idx] = gatherRows<double>(h.cx, rows, H5::PredType::NATIVE_DOUBLE);
            m_my[file_idx] = gatherRows<double>(h.my, rows, H5::PredType::NATIVE_DOUBLE);
            m_cy[file_idx] = gatherRows<double>(h.cy, rows, H5::PredType::NATIVE_DOUBLE);
            m_chi2[file_idx] = gatherRows<double>(h.chi2, rows, H5::PredType::NATIVE_DOUBLE);
            m_ndof[file_idx] = gatherRows<int32_t>(h.ndof, rows, H5::PredType::NATIVE_INT32);
        }

        for (size_t k = 0; k < entries.size(); ++k) {
            const auto& e = entries[k];
            out[k] = TrackFields{m_mx.at(e.file_idx).at(e.row_idx), m_cx.at(e.file_idx).at(e.row_idx),
                                  m_my.at(e.file_idx).at(e.row_idx), m_cy.at(e.file_idx).at(e.row_idx),
                                  m_chi2.at(e.file_idx).at(e.row_idx), m_ndof.at(e.file_idx).at(e.row_idx)};
        }
        return out;
    }

    void Reader::run(DataBatch& batch) {
        auto& geo = GeometryManager::getInstance();

        size_t b_size = batch_size_;
        if (batch.requested_batch_size > 0) b_size = std::min(b_size, batch.requested_batch_size);

        size_t start = track_cursor_;
        size_t end = std::min(start + b_size, track_index_.size());

        // Builds a Cluster from a fetched row - used for both track-linked
        // clusters and replayed "unlinked" DUT clusters below. Adds a
        // single synthetic Pixel at the cluster's centroid: only the
        // aggregate (col, row) survives the save/replay round-trip, not
        // the individual constituent pixels, but DUTAssociation's default
        // cut metric is "distance to the nearest pixel within the
        // cluster" - without at least one pixel present, that loop never
        // runs and the distance silently stays at its DBL_MAX sentinel,
        // failing every cut unconditionally regardless of true proximity.
        // A centroid stand-in loses the (usually small) nearest-pixel-vs-
        // centroid distinction but keeps the cut meaningful.
        auto makeCluster = [&](const ClusterRow& cr, const std::string& det_name) {
            auto cluster = std::make_shared<Cluster>();
            cluster->setEventID(cr.event_id);
            cluster->setLocalCenter(cr.col, cr.row);
            cluster->setGlobal({cr.global_x, cr.global_y, cr.global_z});
            cluster->setTimestamp(cr.timestamp);
            if (geo.hasDetector(det_name)) cluster->setDetector(&geo.getDetector(det_name));
            auto px = std::make_shared<Pixel>(det_name, static_cast<int>(std::lround(cr.col)),
                                               static_cast<int>(std::lround(cr.row)), cr.charge, cr.timestamp,
                                               cr.event_id, cr.channel);
            cluster->addPixel(px);
            return cluster;
        };

        auto detNameOf = [&](uint64_t cluster_id) -> std::string {
            if (cluster_id >= cluster_index_.size() || !cluster_index_[cluster_id].valid) return "";
            return cluster_det_names_[cluster_index_[cluster_id].det_idx];
        };

        uint64_t last_event = 0;
        if (start < end) {
            last_event = track_index_[end - 1].event_id;
            // Delta from the running high-water mark (#last_reported_event_),
            // not this batch's own internal span - see the member's docs
            // for why: it's what lets an untracked-event gap that straddles
            // a batch boundary still get counted instead of silently lost.
            // Accumulated (+=), not assigned - a pipeline can carry more
            // than one Reader instance (one [[Reader]] block per run, see
            // Reader.hpp's class docs), each writing into this same shared
            // DataBatch within one iteration; assignment would let the
            // last instance to run silently clobber any earlier instance's
            // contribution to this batch's event count.
            auto last_event_i = static_cast<int64_t>(last_event);
            batch.events_in_batch += static_cast<size_t>(last_event_i - last_reported_event_);
            last_reported_event_ = last_event_i;

            // Batch-fetch this window's fit parameters and every cluster
            // its /track_cluster_links entries reference in two combined
            // reads (rather than one HDF5 call per row) - see
            // fetchTrackFields()/fetchClusters() for how the scattered
            // rows get coalesced into a handful of hyperslab reads.
            std::vector<TrackIndex> window(track_index_.begin() + static_cast<long>(start),
                                            track_index_.begin() + static_cast<long>(end));
            auto fields = fetchTrackFields(window);

            std::vector<uint64_t> needed_cluster_ids;
            for (const auto& tr : window) {
                auto [ids, count] = clustersForTrack(tr.track_id);
                needed_cluster_ids.insert(needed_cluster_ids.end(), ids, ids + count);
            }
            auto cluster_rows = fetchClusters(needed_cluster_ids);

            for (size_t i = 0; i < window.size(); ++i) {
                const auto& tr = window[i];
                const auto& f = fields[i];
                auto track = std::make_shared<Track>();
                track->setEventID(tr.event_id);
                track->setParams({f.mx, f.cx, f.my, f.cy});
                track->setChi2(f.chi2);
                track->setNdof(f.ndof);

                auto [link_ids, link_count] = clustersForTrack(tr.track_id);
                for (size_t k = 0; k < link_count; ++k) {
                    uint64_t cluster_id = link_ids[k];
                    auto cr_it = cluster_rows.find(cluster_id);
                    if (cr_it == cluster_rows.end()) continue;
                    std::string det_name = detNameOf(cluster_id);
                    if (det_name.empty()) continue;
                    auto cluster = makeCluster(cr_it->second, det_name);
                    track->addCluster(cluster);
                    batch.clusters[det_name].push_back(cluster);
                }

                batch.tracks.push_back(track);
            }
            track_cursor_ = end;
        }

        if (track_index_.empty()) {
            // Pure waveform/cluster replay (no tracks file loaded) - there's
            // no track-batch boundary to synchronize against, so chunk by
            // an advancing SYNTHETIC event-id window instead (b_size events
            // per call), same convention as the tracks-based path below.
            //
            // Loops internally (like EventLoaderHDF5::run()'s own file-
            // boundary fix) rather than returning an empty batch just
            // because THIS particular window happened to contain no rows
            // for any detector - an empty batch is FrameworkManager's own
            // signal for "no more data anywhere", which would wrongly halt
            // the run early on a sparse event-id gap.
            uint64_t start_event = static_cast<uint64_t>(last_reported_event_ + 1);
            while (start_event <= max_loaded_event_id_) {
                uint64_t end_event = std::min(start_event + b_size, max_loaded_event_id_ + 1);
                bool emitted_anything = false;

                for (auto& [det_name, idx] : waveforms_by_det_) {
                    size_t& cursor = waveform_cursor_[det_name];
                    if (cursor >= idx.size() || idx[cursor].event_id >= end_event) continue;

                    size_t window_end = cursor;
                    while (window_end < idx.size() && idx[window_end].event_id < end_event) ++window_end;

                    std::vector<WaveformIndex> wwindow(idx.begin() + static_cast<long>(cursor),
                                                         idx.begin() + static_cast<long>(window_end));
                    auto rows = fetchWaveforms(det_name, wwindow);

                    auto& wf_vec = batch.waveforms[det_name];
                    for (const auto& r : rows) {
                        auto wf = std::make_shared<Waveform>(det_name, r.col, r.row, r.timestamp,
                                                              std::vector<double>{}, r.event_id, r.channel);
                        wf->baseline = r.baseline; wf->noise = r.noise; wf->amplitude = r.amplitude;
                        wf->snr = r.snr; wf->rise_time = r.rise_time;
                        for (size_t k = 0; k < r.start_times.size() && k < cfd_fractions_.size(); ++k) {
                            wf->start_time_by_cfd.emplace_back(cfd_fractions_[k], r.start_times[k]);
                        }
                        wf->start_time = r.start_times.empty() ? 0.0 : r.start_times[0];
                        wf->integral = r.integral; wf->t0 = r.t0; wf->dt = r.dt;
                        wf->time_over_threshold = r.time_over_threshold; wf->time_over_noise = r.time_over_noise;
                        wf->is_processed = true;
                        wf_vec.push_back(wf);
                        emitted_anything = true;
                    }
                    cursor = window_end;
                }
                for (auto& [det_name, idx] : unlinked_clusters_by_det_) {
                    size_t& cursor = cluster_cursor_[det_name];
                    if (cursor >= idx.size() || idx[cursor].event_id >= end_event) continue;

                    size_t window_end = cursor;
                    while (window_end < idx.size() && idx[window_end].event_id < end_event) ++window_end;

                    std::vector<uint64_t> cids;
                    for (size_t k = cursor; k < window_end; ++k) cids.push_back(idx[k].cluster_id);
                    auto cluster_rows = fetchClusters(cids);

                    auto& cl_vec = batch.clusters[det_name];
                    for (uint64_t cid : cids) {
                        auto it = cluster_rows.find(cid);
                        if (it == cluster_rows.end()) continue;
                        cl_vec.push_back(makeCluster(it->second, det_name));
                        emitted_anything = true;
                    }
                    cursor = window_end;
                }

                // Accumulated (+=), not assigned - see the tracks-based
                // path above for why.
                batch.events_in_batch += end_event - start_event;
                last_reported_event_ = static_cast<int64_t>(end_event - 1);

                if (emitted_anything) break;
                start_event = end_event; // this window was empty - try the next one
            }
        } else {
            // Emit only the waveform/cluster rows whose event falls within
            // THIS batch's track event range, walking each detector's own
            // cursor forward - dumping everything in one shot would starve
            // every batch after the first of any DUT data, since only
            // batch 1's tracks ever got a chance to see it. Once tracks are
            // exhausted, flush whatever's left regardless of event id, so
            // trailing waveform/cluster-only events aren't dropped.
            bool tracks_exhausted = (track_cursor_ >= track_index_.size());
            for (auto& [det_name, idx] : waveforms_by_det_) {
                size_t& cursor = waveform_cursor_[det_name];
                // Only touch batch.waveforms[det_name] when there's actually
                // a row to emit - operator[] default-constructs an empty
                // vector on mere access, and a batch full of empty vectors
                // still counts as "not empty" for FrameworkManager's
                // end-of-data check, which would otherwise loop forever
                // once every detector's data is exhausted.
                if (cursor >= idx.size() || (!tracks_exhausted && idx[cursor].event_id > last_event)) continue;

                size_t window_end = cursor;
                while (window_end < idx.size() && (tracks_exhausted || idx[window_end].event_id <= last_event)) {
                    ++window_end;
                }

                std::vector<WaveformIndex> wwindow(idx.begin() + static_cast<long>(cursor),
                                                     idx.begin() + static_cast<long>(window_end));
                auto rows = fetchWaveforms(det_name, wwindow);

                auto& wf_vec = batch.waveforms[det_name];
                for (const auto& r : rows) {
                    auto wf = std::make_shared<Waveform>(det_name, r.col, r.row, r.timestamp, std::vector<double>{},
                                                          r.event_id, r.channel);
                    wf->baseline = r.baseline; wf->noise = r.noise; wf->amplitude = r.amplitude;
                    wf->snr = r.snr; wf->rise_time = r.rise_time;
                    for (size_t k = 0; k < r.start_times.size() && k < cfd_fractions_.size(); ++k) {
                        wf->start_time_by_cfd.emplace_back(cfd_fractions_[k], r.start_times[k]);
                    }
                    wf->start_time = r.start_times.empty() ? 0.0 : r.start_times[0];
                    wf->integral = r.integral; wf->t0 = r.t0; wf->dt = r.dt;
                        wf->time_over_threshold = r.time_over_threshold; wf->time_over_noise = r.time_over_noise;
                    wf->is_processed = true;
                    wf_vec.push_back(wf);
                }
                cursor = window_end;
            }
            for (auto& [det_name, idx] : unlinked_clusters_by_det_) {
                size_t& cursor = cluster_cursor_[det_name];
                if (cursor >= idx.size() || (!tracks_exhausted && idx[cursor].event_id > last_event)) continue;

                size_t window_end = cursor;
                while (window_end < idx.size() && (tracks_exhausted || idx[window_end].event_id <= last_event)) {
                    ++window_end;
                }

                std::vector<uint64_t> cids;
                for (size_t k = cursor; k < window_end; ++k) cids.push_back(idx[k].cluster_id);
                auto cluster_rows = fetchClusters(cids);

                auto& cl_vec = batch.clusters[det_name];
                for (uint64_t cid : cids) {
                    auto it = cluster_rows.find(cid);
                    if (it == cluster_rows.end()) continue;
                    cl_vec.push_back(makeCluster(it->second, det_name));
                }
                cursor = window_end;
            }
        }
    }

    double Reader::getProgress() const {
        if (track_index_.empty()) {
            // No tracks to size progress against (see run()'s own
            // event-id-windowed chunking for this case) - use the running
            // high-water-mark event id against the highest one loaded
            // instead. max_loaded_event_id_ == 0 with nothing loaded at
            // all is indistinguishable from "just the single event 0" -
            // reads as done either way, correct for both.
            if (max_loaded_event_id_ == 0) return 1.0;
            return std::min(1.0, static_cast<double>(last_reported_event_ + 1) /
                                      static_cast<double>(max_loaded_event_id_ + 1));
        }
        return static_cast<double>(track_cursor_) / static_cast<double>(track_index_.size());
    }

    void Reader::finalize() {
        WR_LOG(STATUS, "Reader finished. Emitted " + std::to_string(track_cursor_) + "/" +
                    std::to_string(track_index_.size()) + " tracks.");
    }

    REGISTER_MODULE(Reader)
}
