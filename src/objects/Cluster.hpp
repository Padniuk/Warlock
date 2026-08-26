#ifndef WARLOCK_CLUSTER_HPP
#define WARLOCK_CLUSTER_HPP

#include "Pixel.hpp"
#include "core/detector/GeometryManager.hpp"
#include <Eigen/Dense>
#include <vector>
#include <memory>
#include <string>
#include <algorithm>

namespace framework {
    class Cluster : public Object {
    public:
        Cluster() : Object(0) {}

        void addPixel(std::shared_ptr<Pixel> p) { 
            pixels_.push_back(p); 
            charge_ += p->charge(); 
        }
        
        const std::vector<std::shared_ptr<Pixel>>& pixels() const { return pixels_; }
        
        void setLocalCenter(double c, double r) { col_ = c; row_ = r; }
        double column() const { return col_; }
        double row() const { return row_; }

        void setGlobal(Eigen::Vector3d p) { pos_ = std::move(p); }
        const Eigen::Vector3d& global() const { return pos_; }

        void setDetector(DetectorGeo* det) {
            detector_ = det; 
            det_id_ = det->name; 
        }
        DetectorGeo* detector() const { return detector_; }
        const std::string& detectorID() const { return det_id_; }

        double charge() const { return charge_; }
        size_t size() const { return pixels_.size(); }

        // Number of distinct columns/rows spanned by this cluster's pixels
        // (matches corry's Cluster::columnWidth()/rowWidth()).
        int columnWidth() const {
            if (pixels_.empty()) return 0;
            int lo = pixels_[0]->column(), hi = lo;
            for (auto const& p : pixels_) { lo = std::min(lo, p->column()); hi = std::max(hi, p->column()); }
            return hi - lo + 1;
        }
        int rowWidth() const {
            if (pixels_.empty()) return 0;
            int lo = pixels_[0]->row(), hi = lo;
            for (auto const& p : pixels_) { lo = std::min(lo, p->row()); hi = std::max(hi, p->row()); }
            return hi - lo + 1;
        }

        // Timestamp of the earliest-arriving pixel in the cluster (matches
        // corry's Cluster::getEarliestPixel()->timestamp()), distinct from
        // timestamp()/Object's own charge-weighted-average timestamp set by
        // ClusteringSpatial. Identical to that average for a single-pixel
        // cluster; only differs for a multi-pixel cluster.
        double earliestPixelTimestamp() const {
            if (pixels_.empty()) return timestamp();
            double t = pixels_[0]->timestamp();
            for (auto const& p : pixels_) t = std::min(t, p->timestamp());
            return t;
        }

        // Raw hardware channel (Pixel::channel()) of this cluster's
        // highest-charge pixel, matching the "seed pixel" convention used
        // elsewhere (see ClusteringSpatial's clusterSeedCharge). -1 if the
        // cluster has no pixels, or if its pixel never had real channel
        // info (e.g. reconstructed from a saved file).
        int seedChannel() const {
            if (pixels_.empty()) return -1;
            auto seed = pixels_[0];
            for (auto const& p : pixels_) if (p->charge() > seed->charge()) seed = p;
            return seed->channel();
        }

    private:
        DetectorGeo* detector_{nullptr};
        std::string det_id_;
        std::vector<std::shared_ptr<Pixel>> pixels_;
        Eigen::Vector3d pos_{0,0,0};
        double charge_{0};
        double col_{0}, row_{0};
    };

    /// Maps a raw CAEN DT5742 hardware channel to its trigger group - a
    /// pure function of the channel number, per EUDAQ2's own CAEN plugin
    /// (eudaq/user/caen-dt5742/module/src/CAENDT5742RawEvent2StdEventConverter.cc):
    /// channels 0-7 are group 0, 8-15 are group 1, channel 16 is group 0's
    /// own trigger-timing channel ("TR0"), channel 17 is group 1's.
    /// Comparing timing across channels wired to different trigger groups
    /// isn't physically meaningful, since each group samples independently
    /// - relevant when comparing specific channels/clusters across
    /// detectors, e.g. in AnalysisTiming. Returns -1 for channel < 0 or > 17.
    inline int triggerGroupOf(int channel) {
        if (channel < 0 || channel > 17) return -1;
        if (channel < 8) return 0;
        if (channel < 16) return 1;
        return channel == 16 ? 0 : 1;
    }
}
#endif