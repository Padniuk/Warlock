/**
 * @file ClusteringSpatial.hpp
 * @brief Groups adjacent pixel hits into clusters, per detector and per event.
 */

#ifndef WARLOCK_CLUSTERINGSPATIAL_HPP
#define WARLOCK_CLUSTERINGSPATIAL_HPP

#include "core/Module.hpp"
#include <mutex>
#include <string>
#include <map>
#include <set>
#include <vector>

namespace framework {
    struct DetectorGeo;
    /**
     * @brief Flood-fills 8-connected neighboring pixels within each event
     * into Cluster objects (charge-weighted or geometric centroid, per
     * `charge_weighting`), computing each cluster's local/global position
     * from the result. Mirrors corryvreckan's ClusteringSpatial module.
     * Parallelized per detector, further chunked by event range within
     * each detector so the configured thread count is actually used
     * regardless of how many detectors are present.
     */
    class ClusteringSpatial : public Module {
    public:
        ClusteringSpatial(Configuration cfg, Configuration g_cfg, ThreadPool* pool);

        void initialize() override;
        void run(DataBatch& batch) override;
        void finalize() override;

    private:
        bool use_weighting_; ///< Whether cluster centroids are charge-weighted (vs. plain geometric mean).
        /// If non-empty, an event is only clustered if every listed
        /// detector has at least one pixel in it (intersection of their
        /// event sets) - otherwise every detector's own events are used.
        std::vector<std::string> discard_events_missing_on_;
        /// Detector names and/or types left out of clustering entirely -
        /// same exclude-list convention as Tracking4D::exclude_detectors_/
        /// AlignmentMillepede::exclude_detectors_/EventLoaderHDF5::
        /// exclude_detectors_; see #isExcluded().
        std::vector<std::string> exclude_detectors_;
        /// True if `det` matches #exclude_detectors_ by exact name or by
        /// type - mirrors Tracking4D::isExcluded().
        bool isExcluded(const DetectorGeo& det) const;
        std::mutex cluster_mtx_; ///< Guards merging each event-chunk's clusters into the shared DataBatch.
    };
}
#endif
