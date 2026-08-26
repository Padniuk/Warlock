#ifndef WARLOCK_EVENTLOADERHDF5_HPP
#define WARLOCK_EVENTLOADERHDF5_HPP

#include "core/Module.hpp"
#include <H5Cpp.h>
#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <set>
#include <mutex>

namespace framework {
    struct DetectorGeo;

    /**
     * @brief Loads pixel hits and (optionally) waveform samples from an
     * HDF5 file produced by minionRawConverter, the event source at the
     * start of every Warlock pipeline.
     *
     * Reads the file's `/hits` dataset(s), builds a route table mapping
     * each raw sensor identity to the matching detector(s) in the current
     * geometry (a raw sensor can fan out to several Warlock detectors via
     * DetectorGeo::source_sensor / DetectorGeo::source_channels, e.g. one
     * physical multi-die CAEN board split into independent detectors), and
     * emits DataBatch::pixels (always) and DataBatch::waveforms (if
     * `load_waveforms` is set and the file has per-hit sample data)
     * batch by batch. `type`/`exclude_detectors` restrict which detectors
     * get loaded at all, same convention as MaskCreator/Prealignment/
     * AlignmentTrackChi2/DUTAssociation.
     */
    class EventLoaderHDF5 : public Module {
    public:
        EventLoaderHDF5(Configuration cfg, Configuration global_cfg, ThreadPool* pool);
        
        void initialize() override;
        void run(DataBatch& batch) override;
        void finalize() override;

        double getProgress() const override {
            if (file_names_.empty()) return 0.0;
            double file_prog = (total_hits_ > 0) ? (static_cast<double>(current_hit_idx_) / total_hits_) : 0.0;
            return (static_cast<double>(current_file_idx_) + file_prog) / file_names_.size();
        }

    private:
        bool openNextFile();

        std::vector<std::string> file_names_;
        std::vector<size_t> max_events_;
        size_t current_file_idx_{0};

        /// Restricts buildRouteTable() to detectors of this type ("" = every
        /// detector in the .geo file) - same single-string `type`
        /// convention MaskCreator/Prealignment/AlignmentTrackChi2/
        /// DUTAssociation use, replicating corry's own per-module `type`
        /// filter (ModuleManager.cpp).
        std::string type_filter_;
        /// Detector names and/or types left out of loading, independent of
        /// #type_filter_ - same exclude-list convention as
        /// %Tracking4D::exclude_detectors_ / %AlignmentMillepede::exclude_detectors_. Needed
        /// because #type_filter_ can only express a single included type,
        /// but a stage may need several types loaded together minus one
        /// (e.g. "mimosa26" + "cmsit", not "caendt5742").
        std::vector<std::string> exclude_detectors_;
        /// True if `det` matches #exclude_detectors_ by exact name or by
        /// type - mirrors Tracking4D::isExcluded()/AlignmentMillepede::
        /// isExcluded().
        bool isExcluded(const DetectorGeo& det) const;

        size_t current_hit_idx_{0};
        size_t total_hits_{0};

        /// Cumulative event_id offset applied to the file CURRENTLY open in
        /// #h_event_ids_ (0 for the first file). Each raw-converted file's
        /// own /hits/event_id numbering independently starts near 0 (the
        /// DAQ resets its counter per run), so when `file_name` lists
        /// several files as one combined dataset, their ranges would
        /// otherwise collide - a saved, reloaded (Reader) multi-file
        /// output sorts globally by event_id, so a collision would
        /// cross-associate tracks/clusters/waveforms from different files.
        /// Every hit's event_id is shifted by this offset as soon as it's
        /// read (see openNextFile()), making event_id strictly increasing
        /// and unique across the whole multi-file run. Kept here (not just
        /// applied and forgotten) so run()'s per-file number_of_events cap
        /// can be translated into this file's own portion of the global
        /// numbering.
        uint64_t current_file_event_offset_{0};
        /// Offset to apply to the NEXT file opened - #current_file_event_offset_
        /// + 1 past the highest event_id seen in the file just finished.
        uint64_t next_file_event_offset_{0};
        
        H5::H5File current_file_;
        std::unordered_map<uint32_t, std::string> sensor_map_;
        /// Raw CAEN channel -> DC-offset baseline (Volts), from the "/hits"
        /// group's "caen_dc_offset_volts_ch<N>" attributes (written by
        /// minionRawConverter). Attached to each constructed Waveform (see
        /// Waveform::dc_offset_volts) so WaveformProcessingCAEN can decode
        /// "/hits/samples"' stored ADC code back to Volts without reading
        /// HDF5 attributes itself.
        std::unordered_map<int, double> channel_dc_offset_volts_;

        /// One Warlock detector that consumes a given raw identity's hits -
        /// an empty `filter_channels` means "every hit" (an ordinary,
        /// non-split detector, the historical/default case); otherwise
        /// only hits whose raw channel (see "/hits/channel") is listed in
        /// `filter_channels` belong to this detector, and get remapped to
        /// row 0 in its own local frame (see
        /// DetectorGeo::source_sensor/source_channels).
        struct DetectorRoute {
            std::string name;
            std::vector<int> filter_channels;
        };
        /// Raw identity ("type_planeID", e.g. "CAEN_UZH_3") -> every
        /// Warlock detector fed by it - usually one entry (the detector of
        /// the same name), but more than one for a raw sensor split across
        /// several Warlock detectors by channel. Built once, lazily, from
        /// every DetectorGeo's source_sensor/source_channels (or its own
        /// name, for a non-split detector).
        std::unordered_map<std::string, std::vector<DetectorRoute>> route_table_;
        void buildRouteTable();

        /// Raw identities present in the .geo file but excluded by
        /// #type_filter_ - lets the per-hit "unknown sensor" warning in
        /// run() distinguish "filtered by type" (expected, quiet) from
        /// "not in your .geo file at all" (an actual config problem).
        std::set<std::string> filtered_by_type_;

        // Resolved "type_planeID" -> the route_table_ entry to use (nullptr =
        // unknown/not in geometry), keyed by (sensor_id << 32 | plane_id).
        // Avoids rebuilding the name string and re-querying route_table_ for
        // every single hit - both are otherwise redundant work since a
        // handful of (sensor_id, plane_id) combinations recur across every
        // hit. The actual per-hit detector still depends on that hit's own
        // row when the matched entry has more than one route (a split
        // sensor) - see run().
        std::unordered_map<uint64_t, const std::vector<DetectorRoute>*> route_cache_;

        std::vector<uint64_t> h_event_ids_;
        std::vector<int32_t> h_cols_, h_rows_;
        /// Raw hardware channel per hit, -1 for any sensor that doesn't
        /// carry one (every non-CAEN sensor) - see DetectorGeo::source_channels.
        /// Optional: older files predating this column read as all -1.
        std::vector<int32_t> h_channels_;
        std::vector<uint32_t> h_sensor_ids_, h_plane_ids_;
        std::vector<double> h_timestamps_, h_wf_t0_, h_wf_dt_;

        H5::DataSet ds_samples_;
        bool has_waveforms_{false};
        bool load_waveforms_{false};

        // minionRawConverter output stores waveform data compactly: one row
        // per hit that actually has a waveform (see has_waveform), not one
        // padded row per hit overall - most detectors in a mixed
        // telescope+DUT run never produce waveforms, so a 1:1 layout would
        // waste the bulk of the file on zeroed rows. compact_waveform_format_
        // is set once /hits/has_waveform is found to exist; when it doesn't
        // (older files), everything falls back to the full-length, 1:1
        // index-aligned h_wf_t0_/h_wf_dt_/ds_samples_ read below.
        bool compact_waveform_format_{false};
        std::vector<uint8_t> h_has_waveform_;
        std::vector<double> h_wf_t0_compact_, h_wf_dt_compact_;
        size_t waveform_rows_read_{0};
        /// Sensor types already warned about in the legacy (non-compact)
        /// waveform-association fallback - printed once per type, not
        /// once per hit, since that check runs in a per-hit hot loop.
        std::set<std::string> warned_unknown_wf_types_;

        std::set<std::string> unknown_names_;
        std::mutex batch_mtx_;
    };
}
#endif