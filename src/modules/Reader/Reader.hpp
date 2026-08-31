/**
 * @file Reader.hpp
 * @brief Reloads tracks/clusters/waveforms saved by Tracking4D or WaveformProcessingCAEN back into a DataBatch.
 */

#ifndef WARLOCK_READER_HPP
#define WARLOCK_READER_HPP

#include "core/Module.hpp"
#include <H5Cpp.h>
#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace framework {
    /**
     * @brief Loads fitted tracks/clusters (as written by Tracking4D's
     * optional storage_file) and per-channel waveform parameters (as
     * written by WaveformProcessingCAEN's optional storage_file) back into a
     * DataBatch, so downstream modules (WaveformSelector, DUTAssociation,
     * ...) can reuse already-fitted data instead of re-running
     * clustering/tracking from raw hits every time. Takes one config
     * parameter, `filenames` (an array - a single file works fine as a
     * one-element array, see Configuration::getArray()), and reads
     * whichever of /tracks, /track_cluster_links, /clusters/\<name\>,
     * /waveforms/\<name\> each file actually contains.
     *
     * Combining several independent runs: don't list two runs' files
     * together in one Reader's `filenames` - event_id/track_id/cluster_id
     * are each assigned independently per run (typically starting near 0
     * every time, same as raw hit event_id before EventLoaderHDF5 offsets
     * it), so two different runs' values collide routinely and a single
     * Reader instance has no way to tell them apart. Instead, add one
     * `[[Reader]]` block per run (the same array-of-tables convention
     * `[[WaveformSelector]]` uses for one instance per detector - see
     * FrameworkManager::loadConfiguration()) - each instance gets its own completely
     * independent set of member state below, so there is no shared
     * structure for two runs' ids to collide in. Every instance's `run()`
     * writes into the same shared DataBatch each pipeline iteration, and
     * the framework's end-of-data check already ANDs emptiness across
     * every module, so this composes correctly with no Reader-side
     * multi-run logic needed at all.
     *
     * Memory model: only cheap per-row *index* columns (event_id, track_id,
     * cluster_id - a few bytes each) are read fully into memory in
     * initialize(). The heavy per-row payload (track fit parameters,
     * cluster coordinates, waveform samples) is left on disk - source HDF5
     * datasets are kept open for the module's lifetime and fetched via
     * targeted, batched hyperslab reads scoped to whatever window of rows
     * the current run() call actually needs, then discarded. This bounds
     * peak memory to roughly one batch's worth of hydrated rows (matching
     * EventLoaderHDF5's own "batch_size events resident, then freed"
     * model) instead of the entire file's content, regardless of on-disk
     * row ordering (see the .cpp for why on-disk order can't be assumed
     * sorted, and how the index-then-gather design sidesteps that).
     */
    class Reader : public Module {
    public:
        Reader(Configuration cfg, Configuration g_cfg, ThreadPool* pool);

        void initialize() override;
        void run(DataBatch& batch) override;
        void finalize() override;
        double getProgress() const override;

    private:
        // ---- Lightweight per-row index entries, kept fully in memory ----
        // (a few bytes/row - cheap even at tens of millions of rows; the
        // heavy fields these point at are fetched on demand, see below).
        struct TrackIndex {
            uint64_t event_id, track_id;
            uint32_t file_idx, row_idx; ///< Where to fetch this row's heavy fields from.
        };
        struct WaveformIndex {
            uint64_t event_id;
            uint32_t file_idx, row_idx;
        };
        /// Where cluster_id's heavy fields live, plus its event_id (cheap
        /// to carry here since it has to be read to build this index
        /// anyway, and doing so up front means #unlinked_clusters_by_det_
        /// never needs a separate disk read just to sort by event_id).
        /// cluster_id is a dense, per-file-monotonic counter (see .cpp),
        /// so this is indexed directly by cluster_id rather than kept in a
        /// map - at tens of millions of clusters this is itself a
        /// meaningful chunk of memory, so field order here is deliberate
        /// (event_id first, so the struct's own required 8-byte alignment
        /// doesn't leave padding gaps between smaller fields): 24 bytes/
        /// entry this way instead of 32 with the naive declaration order.
        struct ClusterLoc {
            uint64_t event_id;
            uint32_t file_idx, row_idx;
            uint16_t det_idx;
            bool valid{false};
            /// True once this cluster_id has been seen in a
            /// /track_cluster_links row - initialize() uses this instead
            /// of a separate std::set<uint64_t> of every linked id (which
            /// at tens of millions of links was the single largest
            /// allocation in this whole module) to decide which clusters
            /// belong in #unlinked_clusters_by_det_. Costs nothing extra:
            /// shares the padding byte(s) already needed to round this
            /// struct up to its 8-byte alignment.
            bool linked{false};
        };
        /// One entry per cluster never referenced by any /track_cluster_links
        /// row (e.g. RD53B_115's own clusters, saved by
        /// Tracking4D::saveDutClusters()) - kept sorted by event_id per
        /// detector so run() can replay them event-synced, fetching full
        /// rows through the same fetchClusters() path as track-linked ones.
        struct ClusterIndexEntry {
            uint64_t event_id, cluster_id;
        };

        /// A fully-hydrated cluster row, materialized only for whichever
        /// cluster_ids the current batch actually needs.
        struct ClusterRow {
            uint64_t event_id;
            double col, row, global_x, global_y, global_z, charge, timestamp;
            int32_t channel{-1};
        };
        /// A fully-hydrated waveform row, materialized only for whichever
        /// rows fall in the current replay window.
        struct WaveformRow {
            uint64_t event_id;
            int32_t col, row;
            int32_t channel{-1};
            double timestamp, baseline, noise, amplitude, snr, rise_time, integral, t0, dt;
            /// 0.0 (rather than left uninitialized) for a file saved before
            /// these two existed - see WaveformDetHandles::has_tot_ton.
            double time_over_threshold{0.0}, time_over_noise{0.0};
            std::vector<double> start_times;
        };

        /// Open dataset handles for one file's /tracks heavy fields.
        struct TrackFileHandles {
            H5::DataSet mx, cx, my, cy, chi2, ndof;
        };
        /// Open dataset handles for one (detector, file)'s /clusters/<det>
        /// heavy fields. No event_id handle - already captured per-row in
        /// #ClusterLoc at index-build time.
        struct ClusterDetHandles {
            H5::DataSet col, row, global_x, global_y, global_z, charge, timestamp;
            H5::DataSet channel;
            bool has_channel{false};
        };
        /// Open dataset handles for one (detector, file)'s /waveforms/<det> heavy fields.
        struct WaveformDetHandles {
            H5::DataSet col, row, timestamp, baseline, noise, amplitude, snr, rise_time, integral, t0, dt;
            H5::DataSet channel;
            bool has_channel{false};
            /// Optional: only present in files saved after time_over_threshold/
            /// time_over_noise were added - WaveformRow's copies default to
            /// 0.0 for an older file, same fallback convention as #has_channel.
            H5::DataSet time_over_threshold, time_over_noise;
            bool has_tot_ton{false};
            H5::DataSet start_time_cfd; ///< Rank-2 (n_rows x cfd_width), or a legacy rank-1 "start_time" (see #cfd_rank2).
            hsize_t cfd_width{0};
            bool cfd_rank2{false}; ///< False for a legacy rank-1 "start_time" dataset (cfd_width == 1 either way).
        };

        std::vector<std::string> filenames_; ///< Input HDF5 files, loaded (and concatenated) in initialize().
        size_t batch_size_; ///< Tracks emitted per run() call, capped by DataBatch::requested_batch_size. Also the row-fetch granularity for waveform/cluster replay windows.

        std::vector<std::unique_ptr<H5::H5File>> files_; ///< Kept open for the module's lifetime so heavy fields can be fetched on demand.

        // ---- Tracks ----
        std::vector<TrackIndex> track_index_; ///< Sorted by event_id after initialize().
        std::vector<TrackFileHandles> track_handles_; ///< track_handles_[file_idx], only populated for files that have /tracks.
        size_t track_cursor_{0}; ///< Index into #track_index_ of the next row to emit.

        // ---- Track -> cluster links ----
        // Raw (track_id, cluster_id) pairs, appended across every loadFile()
        // call - built into the CSR-style grouping below exactly once, in
        // initialize(), then cleared (see buildTrackClusterIndex()). A
        // std::map<uint64_t, std::vector<uint64_t>> used to hold this
        // permanently instead: at tens of millions of links, a red-black
        // tree's per-node pointer/allocator overhead plus one heap
        // allocation per distinct track_id's own std::vector made it
        // roughly 2.5x heavier than the flat CSR form below for the same
        // data (measured on a real 42M-link file: ~1.1-1.3GB vs ~440MB).
        std::vector<std::pair<uint64_t, uint64_t>> track_cluster_pairs_raw_;

        /// CSR-style track_id -> cluster_ids grouping, built once by
        /// buildTrackClusterIndex(): tc_cluster_ids_[tc_offsets_[i] ..
        /// tc_offsets_[i+1]) are the cluster_ids linked to tc_track_ids_[i]
        /// (sorted, unique). Looked up via #clustersForTrack() - a binary
        /// search over a flat, contiguous array, with none of a tree's
        /// pointer-chasing or a map<vector<>>'s double indirection.
        std::vector<uint64_t> tc_track_ids_;
        std::vector<uint32_t> tc_offsets_;
        std::vector<uint64_t> tc_cluster_ids_;

        /// Sorts #track_cluster_pairs_raw_ and grinds it into the
        /// #tc_track_ids_/#tc_offsets_/#tc_cluster_ids_ CSR form, marking
        /// every referenced cluster_id's ClusterLoc::linked along the way
        /// (see that field's own docs for why). Clears
        /// #track_cluster_pairs_raw_ afterward - it's only ever needed
        /// once, during this call.
        void buildTrackClusterIndex();

        /// Cluster ids linked to `track_id` (a raw pointer into
        /// #tc_cluster_ids_ plus a count, not an owning copy - valid for
        /// this Reader's lifetime), or {nullptr, 0} if `track_id` has no
        /// linked clusters at all.
        std::pair<const uint64_t*, size_t> clustersForTrack(uint64_t track_id) const;

        // ---- Clusters ----
        std::vector<std::string> cluster_det_names_; ///< det_idx -> name.
        std::map<std::string, size_t> cluster_det_idx_; ///< name -> det_idx, build-time only.
        std::vector<ClusterLoc> cluster_index_; ///< Dense, indexed directly by cluster_id.
        /// (det_idx, file_idx) -> open handles - almost always exactly one
        /// file per detector in practice (Reader is normally given one
        /// tracks file + one waveforms file, never two files both
        /// defining the same detector's clusters), but kept general.
        std::map<std::pair<size_t, size_t>, ClusterDetHandles> cluster_handles_;
        std::map<std::string, std::vector<ClusterIndexEntry>> unlinked_clusters_by_det_; ///< Per-detector, sorted by event_id.
        std::map<std::string, size_t> cluster_cursor_;

        // ---- Waveforms ----
        std::map<std::string, std::vector<WaveformIndex>> waveforms_by_det_; ///< Per-detector, sorted by event_id.
        std::map<std::pair<std::string, size_t>, WaveformDetHandles> waveform_handles_; ///< (det_name, file_idx) -> open handles.
        std::map<std::string, size_t> waveform_cursor_;
        /// CFD fractions the loaded start_time_cfd columns correspond to
        /// (same order as WaveformRow::start_times), read once from the
        /// first "cfd_fractions" attribute encountered - every detector
        /// group in a given stage-60 run shares the same WaveformProcessingCAEN
        /// config, so one shared list is correct for the whole file set.
        std::vector<double> cfd_fractions_;

        /// Highest event_id reported to the framework's progress bar so
        /// far (-1 = nothing reported yet). Each run() call sets
        /// `DataBatch::events_in_batch` to the delta between this and the
        /// current batch's own highest event_id, then advances it - a
        /// running high-water mark (rather than each batch's own
        /// `last - first + 1` span) is what lets an untracked-event gap
        /// straddling a batch boundary still get counted, instead of
        /// falling into neither batch's window and being silently lost.
        /// Also used directly as the next synthetic event-id window's
        /// start by the "#track_index_ is empty" path in run() (there's no
        /// track batch to derive one from), and by getProgress() for that
        /// same case.
        int64_t last_reported_event_{-1};
        /// Highest event_id across every loaded track/waveform/cluster
        /// row (computed once in initialize(), after everything is
        /// sorted) - the upper bound the "#track_index_ is empty" pure-
        /// waveform/cluster replay path (no natural per-batch event
        /// boundary otherwise) chunks its synthetic event-id windows up
        /// to in run(), and getProgress() sizes its fraction against.
        uint64_t max_loaded_event_id_{0};

        /// Reads every group this file actually contains, appending onto
        /// the member containers above - safe to call once per configured
        /// filename. Opens (and keeps open) the heavy-field dataset
        /// handles; only cheap index columns are actually read here.
        void loadFile(const std::string& filename, size_t file_idx);

        /// Fetches full rows for exactly the given cluster_ids (any mix of
        /// detectors/files), via one batched hyperslab gather read per
        /// field per (detector, file) group involved. Missing/unresolved
        /// cluster_ids are simply absent from the returned map.
        std::map<uint64_t, ClusterRow> fetchClusters(const std::vector<uint64_t>& cluster_ids);

        /// Fetches full rows for exactly the given #WaveformIndex entries
        /// of one detector (possibly spanning several files), preserving
        /// the caller's entry order in the returned vector.
        std::vector<WaveformRow> fetchWaveforms(const std::string& det_name, const std::vector<WaveformIndex>& entries);

        /// Fetches full fit-parameter rows for exactly the given
        /// #TrackIndex entries (possibly spanning several files),
        /// preserving the caller's entry order.
        struct TrackFields { double mx, cx, my, cy, chi2; int32_t ndof; };
        std::vector<TrackFields> fetchTrackFields(const std::vector<TrackIndex>& entries);
    };
}
#endif
