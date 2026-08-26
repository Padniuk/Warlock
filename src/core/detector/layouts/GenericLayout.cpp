/** @file GenericLayout.cpp */
#include "core/detector/layouts/GenericLayout.hpp"
#include "core/detector/GeometryManager.hpp"
#include "core/utils/Logger.hpp"
#include "core/utils/PlotManager.hpp"

namespace framework {

    double GenericLayout::effectivePitchX(const DetectorGeo& det) const {
        // Same column-spacing formula as every other layout: a genuine 2D
        // grid never has a die gap, so pitch_x_gap is 0 here in practice,
        // but the formula still degrades correctly if it isn't.
        return (det.n_pixels_x >= 2) ? (det.pitch_x + det.pitch_x_gap / (det.n_pixels_x - 1)) : det.pitch_x;
    }

    double GenericLayout::effectivePitchY(const DetectorGeo& det) const {
        // Same row-spacing formula as every other layout: for a
        // single-row board this is just pitch_y; a genuine 2D grid never
        // has a die gap, so pitch_y_gap is 0 here in practice, but the
        // formula still degrades correctly if it isn't.
        return (det.n_pixels_y >= 2) ? (det.pitch_y + det.pitch_y_gap / (det.n_pixels_y - 1)) : det.pitch_y;
    }

    void GenericLayout::prealign(DetectorGeo& det, const Histogram1D& hist_x, const Histogram1D& hist_y) const {
        // A genuine 2D pixel grid is still a single hard-edged rectangular
        // footprint, so the same combined-plateau edge-detection
        // (Layout::fitPlateauCenter(), inherited) as every other layout
        // applies here.
        for (const std::string axis : {"x", "y"}) {
            const Histogram1D& hist = (axis == "x") ? hist_x : hist_y;
            if (hist.getEntries() < 100) continue;

            double width = (axis == "x") ? det.sizeX() : det.sizeY();
            double pitch = (axis == "x") ? det.pitch_x : det.pitch_y;
            auto fit = fitPlateauCenter(hist, width, pitch);
            double center = fit.center;

            // track->positionAt(det_z) evaluates the fitted line in pure
            // global coordinates - it never subtracts the detector's
            // current position, so `center` (the fitted plateau center)
            // already IS the detector's correct global position, not an
            // offset to add on top of whatever nominal/placeholder value
            // it started at.
            if (axis == "x") det.position(0) = center;
            else det.position(1) = center;

            // Not a Module, so no getName()/WR_LOG - logged under
            // "PrealignmentCAEN" directly (Logger::log()), matching where
            // this alignment step conceptually runs from the user's
            // perspective.
            Logger::log(LogLevel::STATUS, "PrealignmentCAEN",
                        det.name + " " + axis + " Aligned -> Center: " + std::to_string(center) +
                            "mm [derived width: " + std::to_string(fit.derived_width) +
                            "mm, expected from geometry: " + std::to_string(width) + "mm]");
            if (!fit.warning.empty()) {
                Logger::log(LogLevel::WARNING, "PrealignmentCAEN", det.name + " " + axis + ": " + fit.warning);
            }
        }
    }

    const GenericLayout& GenericLayout::instance() {
        static const GenericLayout inst;
        return inst;
    }
}
