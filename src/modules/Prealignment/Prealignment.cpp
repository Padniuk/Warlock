#include "Prealignment.hpp"
#include "core/ModuleFactory.hpp"
#include "core/detector/GeometryManager.hpp"
#include "core/utils/PlotManager.hpp"
#include "core/utils/Logger.hpp"
#include "core/utils/Units.hpp"
#include <algorithm>
#include <cmath>

namespace framework {

    Prealignment::Prealignment(Configuration cfg, Configuration g_cfg, ThreadPool* pool)
        : Module(std::move(cfg), std::move(g_cfg), pool) {
        
        method_               = config.get<std::string>("method", "mean");
        range_abs_            = config.get<double>("range_abs", 10.0);
        damping_factor_       = config.get<double>("damping_factor", 1.0);
        max_correlation_rms_  = config.get<double>("max_correlation_rms", 6.0);
        time_cut_abs_         = config.get<double>("time_cut_abs", 1e9); 
        fixed_planes_         = config.getArray<std::string>("fixed_planes");
        type_filter_          = config.get<std::string>("type", "");
    }

    void Prealignment::initialize() {
        auto& geo = GeometryManager::getInstance();
        std::string ref = geo.getReferenceName();
        WR_LOG(STATUS, "Initializing Prealignment...");

        for (const auto& name : geo.getDetectorNames()) {
            auto& det = geo.getDetector(name);
            auto& ref_det = geo.getDetector(ref);

            // The reference detector still gets its own plots registered
            // and filled below, matching corry (which histograms the
            // reference too, only excluding it from the shift itself in
            // finalize() - correlating it against its own clusters would
            // make a shift meaningless, but the histogram is still valid).
            if (!type_filter_.empty() && det.type != type_filter_) continue;
            if (det.role == "auxiliary") continue;

            auto& pm = PlotManager::getInstance();
            std::string dir = getName() + "/" + name;

            // Binning matched exactly to corryvreckan's Prealignment module.
            pm.registerPlot1D(dir, "correlationX", 1000, -range_abs_, range_abs_);
            pm.registerPlot1D(dir, "correlationY", 1000, -range_abs_, range_abs_);

            pm.registerPlot2D(dir, "correlationX_2D", 100, -range_abs_, range_abs_, 100, -range_abs_, range_abs_);
            pm.registerPlot2D(dir, "correlationY_2D", 100, -range_abs_, range_abs_, 100, -range_abs_, range_abs_);
            
            pm.registerPlot2D(dir, "correlationX_2Dlocal", 
                det.n_pixels_x, -0.5, det.n_pixels_x - 0.5, ref_det.n_pixels_x, -0.5, ref_det.n_pixels_x - 0.5);
            pm.registerPlot2D(dir, "correlationY_2Dlocal", 
                det.n_pixels_y, -0.5, det.n_pixels_y - 0.5, ref_det.n_pixels_y, -0.5, ref_det.n_pixels_y - 0.5);
        }
    }

    void Prealignment::run(DataBatch& batch) {
        std::string ref_name = GeometryManager::getInstance().getReferenceName();
        if (batch.clusters.count(ref_name) == 0) return;

        std::map<uint64_t, std::map<std::string, std::vector<std::shared_ptr<Cluster>>>> events;
        for (auto const& entry : batch.clusters) {
            for (auto const& cl : entry.second) events[cl->eventID()][entry.first].push_back(cl);
        }

        auto& pm = PlotManager::getInstance();
        for (auto const& ev_entry : events) {
            auto const& dets = ev_entry.second;
            if (dets.count(ref_name) == 0) continue;
            auto const& ref_cls = dets.at(ref_name);

            for (auto const& [name, clusters] : dets) {
                std::string dir = getName() + "/" + name;
                if (!pm.hasPlot1D(dir, "correlationX")) continue;

                for (auto const& clus : clusters) {
                    for (auto const& ref_cl : ref_cls) {
                        double timeDifference = ref_cl->timestamp() - clus->timestamp();
                        if (std::abs(timeDifference) >= time_cut_abs_) continue;

                        // Reference minus target, matching corry's own convention.
                        double dx = ref_cl->global().x() - clus->global().x();
                        double dy = ref_cl->global().y() - clus->global().y();

                        pm.fill1D(dir, "correlationX", dx);
                        pm.fill1D(dir, "correlationY", dy);
                        
                        pm.fill2D(dir, "correlationX_2D", clus->global().x(), ref_cl->global().x());
                        pm.fill2D(dir, "correlationY_2D", clus->global().y(), ref_cl->global().y());
                        pm.fill2D(dir, "correlationX_2Dlocal", clus->column(), ref_cl->column());
                        pm.fill2D(dir, "correlationY_2Dlocal", clus->row(), ref_cl->row());
                    }
                }
            }
        }
    }

    void Prealignment::finalize() {
        auto& geo = GeometryManager::getInstance();
        std::string ref = geo.getReferenceName();

        for (const auto& name : geo.getDetectorNames()) {
            std::string dir = getName() + "/" + name;
            if (!PlotManager::getInstance().hasPlot1D(dir, "correlationX")) continue;

            auto& hx = PlotManager::getInstance().getPlot1D(dir, "correlationX");
            auto& hy = PlotManager::getInstance().getPlot1D(dir, "correlationY");

            double rmsX = hx.getStdDev();
            double rmsY = hy.getStdDev();
            if (rmsX > max_correlation_rms_ || rmsY > max_correlation_rms_) {
                // Matches corry: this is a warning, not a reason to skip the
                // shift - the detector is still moved using whatever
                // (unreliable) correlation peak was found. Runs
                // unconditionally, including for the reference detector;
                // only the shift computation/application below is skipped
                // for it.
                WR_LOG(ERROR, "Detector " + name + ": RMS is too wide for prealignment shifts (RMS X = " +
                                   Units::display(rmsX) + ", RMS Y = " + Units::display(rmsY) + ").");
            }

            bool is_fixed = std::find(fixed_planes_.begin(), fixed_planes_.end(), name) != fixed_planes_.end();
            if (name == ref || is_fixed) continue;

            double sx = 0, sy = 0;
            if (method_ == "maximum") {
                sx = hx.getMaxBinCenter();
                sy = hy.getMaxBinCenter();
            } else if (method_ == "gauss_fit") {
                sx = hx.getPrecisePeak();
                sy = hy.getPrecisePeak();
            } else {
                sx = hx.getMean();
                sy = hy.getMean();
            }

            auto& det = geo.getDetector(name);
            double dx_final = sx * damping_factor_;
            double dy_final = sy * damping_factor_;

            det.position(0) += dx_final;
            det.position(1) += dy_final;

            WR_LOG(INFO, "Running detector " + name);
            WR_LOG(INFO, "Using prealignment method: " + method_);
            
            WR_LOG(INFO, "Move in x by = " + Units::display(dx_final) + ", and in y by = " + Units::display(dy_final));
            WR_LOG(INFO, "Detector position after shift in x = " + Units::display(det.position(0)) + ", and in y = " + Units::display(det.position(1)));
        }

        std::string out_geo = global_config.get<std::string>("detectors_file_updated", "updated.geo");
        geo.saveGeometry(out_geo);
        WR_LOG(STATUS, "Prealignment finished. New geometry: " + out_geo);
    }

    REGISTER_MODULE(Prealignment)
}