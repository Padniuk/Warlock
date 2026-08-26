/** @file OrthogonalLayout.cpp */
#include "core/detector/layouts/OrthogonalLayout.hpp"
#include "core/detector/GeometryManager.hpp"
#include "core/utils/Logger.hpp"
#include "core/utils/PlotManager.hpp"

namespace framework {

    double OrthogonalLayout::effectivePitchX(const DetectorGeo& det) const {
        // No known orthogonal-split board splits by column, so this is
        // just the generic formula - a future column-split variant would
        // get its own Layout subclass rather than a branch here.
        return (det.n_pixels_x >= 2) ? (det.pitch_x + det.pitch_x_gap / (det.n_pixels_x - 1)) : det.pitch_x;
    }

    double OrthogonalLayout::effectivePitchY(const DetectorGeo& det) const {
        // pitch_y is each die's own single-pad pitch; the two dies sit on
        // the same carrier with an extra physical gap between them
        // (pitch_y_gap, 0 for a board with no real gap e.g. TREF), spread
        // evenly across the n_pixels_y - 1 gaps between consecutive
        // declared rows.
        return (det.n_pixels_y >= 2) ? (det.pitch_y + det.pitch_y_gap / (det.n_pixels_y - 1)) : det.pitch_y;
    }

    void OrthogonalLayout::prealign(DetectorGeo& det, const Histogram1D& hist_x, const Histogram1D& hist_y) const {
        // The two independently-mounted dies still form a single
        // hard-edged rectangular footprint once sizeY() folds in
        // pitch_y_gap, so a combined-plateau edge fit locates the whole
        // board correctly - it's individual pixels/dies that need
        // dieOf(), not this whole-board alignment step.
        for (const std::string axis : {"x", "y"}) {
            const Histogram1D& hist = (axis == "x") ? hist_x : hist_y;
            if (hist.getEntries() < 100) continue;

            double width = (axis == "x") ? det.sizeX() : det.sizeY();
            double pitch = (axis == "x") ? det.pitch_x : det.pitch_y;
            auto fit = fitPlateauCenter(hist, width, pitch);
            double center = fit.center;

            if (axis == "x") det.position(0) = center;
            else det.position(1) = center;

            Logger::log(LogLevel::STATUS, "PrealignmentCAEN",
                        det.name + " " + axis + " Aligned -> Center: " + std::to_string(center) +
                            "mm [derived width: " + std::to_string(fit.derived_width) +
                            "mm, expected from geometry: " + std::to_string(width) + "mm]");
            if (!fit.warning.empty()) {
                Logger::log(LogLevel::WARNING, "PrealignmentCAEN", det.name + " " + axis + ": " + fit.warning);
            }
        }
    }

    const OrthogonalLayout& OrthogonalLayout::instance() {
        static const OrthogonalLayout inst;
        return inst;
    }
}
