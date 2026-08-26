/** @file Layout.cpp */
#include "core/detector/layouts/Layout.hpp"
#include "core/utils/PlotManager.hpp"
#include <algorithm>
#include <cmath>
#include <sstream>

namespace framework {

    std::string Layout::dieOf(const DetectorGeo&, double, double) const {
        return "";
    }

    std::vector<std::string> Layout::dieLabels(const DetectorGeo&) const {
        return {};
    }

    PlateauFit Layout::fitPlateauCenter(const Histogram1D& hist, double width, double pitch,
                                         double width_tolerance) const {
        size_t n = hist.getBins();
        if (n == 0) return {0.0, 0.0, ""};

        size_t peak_bin = 0;
        double peak_val = hist.getBinContent(0);
        for (size_t i = 1; i < n; ++i) {
            double v = hist.getBinContent(i);
            if (v > peak_val) { peak_val = v; peak_bin = i; }
        }

        double bin_width = hist.getBinWidth();

        // A DUT can have multiple pixels along one axis (e.g. a 2x2
        // sensor), each showing up as its own resolved peak with a dip
        // between them rather than one continuous plateau. Stopping at
        // the first threshold crossing out from the tallest peak would
        // land on that dip and undershoot the true extent, so
        // search_margin spans the whole configured width (plus slack for
        // smearing at the true outer edges) and the scan takes the
        // farthest-out bin still above threshold within that window.
        size_t search_margin_bins = static_cast<size_t>(std::round(1.3 * width / bin_width)) + 2;
        size_t bg_margin_bins = static_cast<size_t>(std::round(1.6 * width / bin_width)) + 2;

        auto local_avg = [&](size_t center_bin) {
            size_t lo = (center_bin > 2) ? center_bin - 2 : 0;
            size_t hi = std::min(center_bin + 2, n - 1);
            double sum = 0.0; size_t count = 0;
            for (size_t i = lo; i <= hi; ++i) { sum += hist.getBinContent(i); ++count; }
            return count > 0 ? sum / count : 0.0;
        };

        double plateau_level = local_avg(peak_bin);
        size_t bg_left_bin = (peak_bin > bg_margin_bins) ? peak_bin - bg_margin_bins : 0;
        size_t bg_right_bin = std::min(peak_bin + bg_margin_bins, n - 1);
        double background_level = 0.5 * (local_avg(bg_left_bin) + local_avg(bg_right_bin));
        double threshold = background_level + 0.5 * (plateau_level - background_level);

        size_t search_lo = (peak_bin > search_margin_bins) ? peak_bin - search_margin_bins : 0;
        size_t search_hi = std::min(peak_bin + search_margin_bins, n - 1);

        size_t leftmost = peak_bin;
        for (size_t i = search_lo; i <= peak_bin; ++i) {
            if (hist.getBinContent(i) >= threshold) { leftmost = i; break; }
        }
        size_t rightmost = peak_bin;
        for (size_t i = search_hi; ; --i) {
            if (hist.getBinContent(i) >= threshold) { rightmost = i; break; }
            if (i == peak_bin) break;
        }

        double left_edge = hist.getBinCenter(0);
        if (leftmost > 0) {
            double v_hi = hist.getBinContent(leftmost);
            double v_lo = hist.getBinContent(leftmost - 1);
            double frac = (v_hi != v_lo) ? (threshold - v_lo) / (v_hi - v_lo) : 0.0;
            left_edge = hist.getBinCenter(leftmost - 1) + frac * bin_width;
        } else {
            left_edge = hist.getBinCenter(leftmost);
        }

        double right_edge = hist.getBinCenter(n - 1);
        if (rightmost + 1 < n) {
            double v_hi = hist.getBinContent(rightmost);
            double v_lo = hist.getBinContent(rightmost + 1);
            double frac = (v_hi != v_lo) ? (threshold - v_lo) / (v_hi - v_lo) : 0.0;
            right_edge = hist.getBinCenter(rightmost + 1) - frac * bin_width;
        } else {
            right_edge = hist.getBinCenter(rightmost);
        }

        double center = (left_edge + right_edge) / 2.0;
        double derived_width = right_edge - left_edge;

        // A channel missing entirely from this sample undershoots the
        // configured width substantially; 0.7 comfortably separates that
        // from an intact multi-channel plateau, which typically lands
        // within a few % of the configured width.
        if (width > 0.0 && derived_width < 0.7 * width) {
            // Shift by half a pixel pitch (a fixed, known hardware
            // quantity), not half of (width - derived_width): the latter
            // overstates the true offset when the surviving channel's own
            // fit comes out abnormally narrow.
            double shift = (pitch > 0.0) ? (pitch / 2.0) : ((width - derived_width) / 2.0);

            // A partially-dead channel (cross-talk, noise triggers)
            // usually still leaves some trace beyond the found edge, while
            // a fully dead channel leaves only background scatter, where a
            // single noisy bin could register as "more" on one side by
            // chance. Requiring a contiguous run of several bins above a
            // looser threshold distinguishes a real (even weak) peak from
            // noise.
            double loose_threshold = background_level + 0.15 * (plateau_level - background_level);
            size_t min_run_bins = std::max<size_t>(3, static_cast<size_t>(std::round(0.3 * pitch / bin_width)));

            auto longest_run_excess = [&](size_t lo, size_t hi) -> double {
                double best_excess = 0.0, run_excess = 0.0;
                size_t run_len = 0;
                for (size_t i = lo; i <= hi; ++i) {
                    double v = hist.getBinContent(i);
                    if (v > loose_threshold) {
                        ++run_len;
                        run_excess += (v - background_level);
                    } else {
                        if (run_len >= min_run_bins) best_excess = std::max(best_excess, run_excess);
                        run_len = 0; run_excess = 0.0;
                    }
                }
                if (run_len >= min_run_bins) best_excess = std::max(best_excess, run_excess);
                return best_excess;
            };

            double excess_left = (leftmost > search_lo) ? longest_run_excess(search_lo, leftmost - 1) : 0.0;
            double excess_right = (rightmost + 1 <= search_hi) ? longest_run_excess(rightmost + 1, search_hi) : 0.0;

            std::ostringstream warning;
            warning << "derived width " << derived_width << "mm undershoots the geometry-expected " << width
                    << "mm - a channel looks missing from this sample.";
            if (excess_left > excess_right && excess_left > 0.0) {
                center -= shift;
                warning << " Shifted center by -" << shift
                        << "mm (half a pixel pitch) toward weak residual signal found on the low side.";
            } else if (excess_right > excess_left && excess_right > 0.0) {
                center += shift;
                warning << " Shifted center by +" << shift
                        << "mm (half a pixel pitch) toward weak residual signal found on the high side.";
            } else {
                warning << " No residual signal on either side to tell which one - center NOT corrected, verify this detector/axis manually.";
            }
            return {center, derived_width, warning.str()};
        }

        // General mismatch check, independent of the severe undershoot
        // branch above: derived_width disagreeing with the configured
        // width in either direction (including an overshoot, which the
        // branch above never catches) is surfaced as a warning rather
        // than left to sit unnoticed in a STATUS line.
        if (width > 0.0) {
            double rel_diff = std::abs(derived_width - width) / width;
            if (rel_diff > width_tolerance) {
                std::ostringstream warning;
                warning << "derived width " << derived_width << "mm " << (derived_width > width ? "overshoots" : "undershoots")
                        << " the geometry-expected " << width << "mm by " << (rel_diff * 100.0)
                        << "% (tolerance " << (width_tolerance * 100.0) << "%) - the configured pitch/gap for this "
                           "axis may not match the real data. Center was NOT withheld (still applied), but verify "
                           "this detector/axis's geometry manually.";
                return {center, derived_width, warning.str()};
            }
        }

        return {center, derived_width, ""};
    }
}
