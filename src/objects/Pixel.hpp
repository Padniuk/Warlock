#ifndef WARLOCK_PIXEL_HPP
#define WARLOCK_PIXEL_HPP
#include "Object.hpp"
#include <string>

namespace framework {
    class Pixel : public Object {
    public:
        /// @param det Owning detector's name.
        /// @param c Pixel column.
        /// @param r Pixel row.
        /// @param ch Charge (or, for a synthetic cluster-centroid pixel, the cluster's total charge).
        /// @param ts Timestamp.
        /// @param ev_id Event ID this hit belongs to.
        /// @param channel Raw hardware channel (0-17 for a CAEN DT5742
        /// board), from minionRawConverter's "/hits/channel" column - -1
        /// (default) for any sensor that doesn't carry one, or for a hit
        /// built without that information at all (e.g. Reader replaying a
        /// saved cluster centroid). Lets a module distinguish which
        /// physical die of a multi-pad board a hit/cluster came from
        /// without needing separate detector geometry per die.
        Pixel(std::string det, int c, int r, double ch, double ts, uint64_t ev_id = 0, int channel = -1)
            : Object(ts, ev_id), det_id_(std::move(det)), col_(c), row_(r), charge_(ch), channel_(channel) {}

        const std::string& detectorID() const { return det_id_; }
        int column() const { return col_; }
        int row() const { return row_; }
        double charge() const { return charge_; }
        int channel() const { return channel_; }

    private:
        std::string det_id_;
        int col_, row_;
        double charge_;
        int channel_;
    };
}
#endif