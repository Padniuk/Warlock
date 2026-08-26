/**
 * @file Tracking4D.hpp
 * @brief Straight-line track finding and fitting across pixel detector planes.
 */

#ifndef WARLOCK_TRACKING4D_HPP
#define WARLOCK_TRACKING4D_HPP

#include "core/Module.hpp"
#include "core/utils/AppendableDataset.hpp"
#include <Eigen/Dense>
#include <H5Cpp.h>
#include <memory>
#include <mutex>
#include <string>
#include <vector>
#include <map>

namespace framework {
    struct DetectorGeo;

    /**
     * @brief Seeds candidate tracks from the outermost (in Z) plane pair
     * carrying clusters in an event, extrapolates through intermediate
     * planes to pick up matching clusters within a spatial/time cut, and
     * refits after every match - mirroring corryvreckan's Tracking4D
     * module. Parallelized per event chunk; plot fills happen inside each
     * worker rather than in a later serial pass (see the timing comments
     * in run()). Optionally dumps fitted tracks and their clusters to HDF5
     * for downstream reuse (see #storage_file_).
     */
    class Tracking4D : public Module {
    public:
        Tracking4D(Configuration cfg, Configuration g_cfg, ThreadPool* pool);

        void initialize() override;
        void run(DataBatch& batch) override;
        void finalize() override;

    private:
        double spatial_cut_x_;    ///< Elliptical spatial-match cut, X half-axis.
        double spatial_cut_y_;    ///< Elliptical spatial-match cut, Y half-axis.
        double time_cut_abs_;     ///< Time cut, in ns, applied to both seed pairs and mid-plane matches.
        size_t min_hits_;         ///< Minimum clusters a candidate must carry to be accepted as a track.
        std::vector<std::string> exclude_detectors_; ///< Detector names and/or types left out of tracking; see #isExcluded().
        bool unique_cluster_usage_; ///< If set, tracks sharing a cluster within an event are deduplicated (see run()).

        /// Per-detector minimum pairwise track distance (um) below which
        /// both tracks get marked not-isolated at that plane - matches
        /// corry's `isolated_planes`/`isolation_cut` (paired by index,
        /// same convention as corry's own config parsing). Empty (default)
        /// means isolation is never evaluated, matching corry's own
        /// "isolation_cut_.size() > 0" gate.
        std::map<std::string, double> isolation_cut_;

        /// True if `det` should be left out of tracking - matched against
        /// #exclude_detectors_ by exact detector name OR by detector type,
        /// so a single list can exclude e.g. "RD53B_115" specifically or
        /// every "caendt5742" at once, without touching that detector's
        /// role in the geometry file.
        bool isExcluded(const DetectorGeo& det) const;

        std::mutex track_mtx_;
        std::atomic<size_t> total_tracks_{0};

        /// Optional: dump fitted tracks and the clusters they're built from
        /// to an HDF5 file, so downstream stages (waveform selection, DUT
        /// association, ...) can reload already-fitted tracks instead of
        /// re-running clustering+tracking from raw hits every time. See the
        /// Reader module for the counterpart that loads this back in. Config
        /// key is `storage_file`, matching WaveformProcessingCAEN's own
        /// parameter name for the same "optional HDF5 dump" concept.
        std::string storage_file_;
        bool save_enabled_{false};
        std::unique_ptr<H5::H5File> save_h5_;

        /// Per-detector set of AppendableDatasets backing one
        /// "/clusters/<det_name>" HDF5 group.
        struct ClusterWriter {
            std::unique_ptr<AppendableDataset> event_id, cluster_id, col, row,
                global_x, global_y, global_z, charge, timestamp, channel;
        };
        std::unique_ptr<AppendableDataset> tr_event_id_, tr_track_id_,
            tr_mx_, tr_cx_, tr_my_, tr_cy_, tr_chi2_, tr_ndof_;
        std::unique_ptr<AppendableDataset> link_track_id_, link_cluster_id_;
        std::map<std::string, ClusterWriter> cluster_writers_;
        uint64_t next_track_id_{0};
        uint64_t next_cluster_id_{0};

        /// Lazily creates (or returns the existing) ClusterWriter for `det_name`.
        ClusterWriter& getClusterWriter(const std::string& det_name);
        /// Appends every track in `tracks` (and the clusters they're built
        /// from) to the HDF5 file opened in initialize().
        void saveTracks(const std::vector<std::shared_ptr<Track>>& tracks);
        /// Saves clusters for real pixel detectors excluded from the fit
        /// (e.g. RD53B_115) as plain per-event rows, with no track link -
        /// unlike saveTracks() above, these were never part of any track's
        /// fit, so there's nothing to link them to. DUTAssociation runs
        /// downstream (in the replay/analysis config) for these detectors,
        /// the same way it does for CAEN boards via WaveformSelector's
        /// waveform-derived clusters - keeping "raw data saving" and "DUT
        /// association" as separate stages. Includes waveform-digitizer
        /// (CAEN) detectors too, since their clusters come from
        /// ClusteringSpatial just like any other excluded detector
        /// (WaveformSelector feeds their selected waveforms into
        /// ClusteringSpatial as Pixels upstream of this module).
        void saveDutClusters(const DataBatch& batch);
    };
}
#endif
