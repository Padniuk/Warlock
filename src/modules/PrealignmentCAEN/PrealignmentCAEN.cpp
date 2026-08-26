/** @file PrealignmentCAEN.cpp */
#include "PrealignmentCAEN.hpp"
#include "core/ModuleFactory.hpp"
#include "core/detector/GeometryManager.hpp"
#include "core/utils/PlotManager.hpp"
#include "core/utils/Logger.hpp"
#include "core/utils/Units.hpp"
#include <cmath>
#include <algorithm>
#include <utility>
#include <sstream>
#include <string>

namespace framework {

    PrealignmentCAEN::PrealignmentCAEN(Configuration cfg, Configuration g_cfg, ThreadPool* pool)
        : Module(std::move(cfg), std::move(g_cfg), pool) {
        dut_names_ = config.getArray<std::string>("dut_names");
    }

    void PrealignmentCAEN::initialize() {
        WR_LOG(STATUS, "Initializing PrealignmentCAEN (Edge Profile Analysis)...");
        for (const auto& name : dut_names_) {
            PlotManager::getInstance().registerPlot1D(getName(), "intercept_x_" + name, 500, -10.0, 10.0);
            PlotManager::getInstance().registerPlot1D(getName(), "intercept_y_" + name, 500, -10.0, 10.0);
        }
    }

    void PrealignmentCAEN::run(DataBatch& batch) {
        if (batch.tracks.empty()) return;
        auto& geo = GeometryManager::getInstance();

        for (const auto& name : dut_names_) {
            if (!geo.hasDetector(name)) continue;
            double det_z = geo.getDetector(name).position.z();

            for (const auto& track : batch.tracks) {
                auto intercept = track->positionAt(det_z);
                PlotManager::getInstance().fill1D(getName(), "intercept_x_" + name, intercept.x());
                PlotManager::getInstance().fill1D(getName(), "intercept_y_" + name, intercept.y());
            }
        }
    }

    void PrealignmentCAEN::finalize() {
        auto& geo = GeometryManager::getInstance();

        for (const auto& name : dut_names_) {
            if (!geo.hasDetector(name)) continue;
            auto& det = geo.getDetector(name);

            // How to locate a DUT from its own intercept histograms is a
            // property of its physical footprint shape, so it lives on
            // det's own Layout (see src/core/detector/layouts/) - not
            // every layout needs to fit a single combined plateau the same
            // way TILGAD/TREF do. No layout set is not necessarily an
            // error (a genuinely single-pad/unsplit CAEN board has no need
            // for one), so GenericLayout still runs the same generic
            // edge-detection rather than skipping alignment outright.
            if (det.hasUnrecognizedLayout()) {
                WR_LOG(WARNING, name + ": unrecognized layout '" + det.layout +
                                     "' - using generic edge-detection alignment as a fallback. Add a "
                                     "dedicated Layout subclass (src/core/detector/layouts/) if this layout "
                                     "needs different handling.");
            }
            auto& hist_x = PlotManager::getInstance().getPlot1D(getName(), "intercept_x_" + name);
            auto& hist_y = PlotManager::getInstance().getPlot1D(getName(), "intercept_y_" + name);
            det.getLayout()->prealign(det, hist_x, hist_y);
        }

        std::string out_geo = global_config.get<std::string>("detectors_file_updated", "caen_aligned.geo");
        geo.saveGeometry(out_geo);
        WR_LOG(STATUS, "PrealignmentCAEN complete. Updated geometry: " + out_geo);
    }

    REGISTER_MODULE(PrealignmentCAEN)
}