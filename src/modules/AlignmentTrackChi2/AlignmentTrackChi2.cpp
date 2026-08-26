/** @file AlignmentTrackChi2.cpp */
#include "AlignmentTrackChi2.hpp"
#include "core/ModuleFactory.hpp"
#include "core/detector/GeometryManager.hpp"
#include "core/utils/PlotManager.hpp"
#include "core/utils/Logger.hpp"
#include "core/utils/Units.hpp"
#include <future>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <unordered_set>

#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>

namespace framework {

    /**
     * @brief Eigen NumericalDiff-compatible functor for one detector's LM
     * fit: `operator()` applies a candidate position/orientation offset,
     * re-projects that detector's clusters, refits every track touching it,
     * and returns each track's sqrt(chi2) as one residual. active_dofs maps
     * the LM parameter vector's indices onto position/orientation
     * components (0-2 = X/Y/Z position, 3-5 = X/Y/Z rotation), so only the
     * DOFs enabled via AlignmentTrackChi2::align_position_ /
     * AlignmentTrackChi2::align_orientation_ are actually varied.
     */
    struct AlignmentFunctor {
        using Scalar = double;
        enum { InputsAtCompileTime = Eigen::Dynamic, ValuesAtCompileTime = Eigen::Dynamic };
        using InputType = Eigen::VectorXd;
        using ValueType = Eigen::VectorXd;
        using JacobianType = Eigen::MatrixXd;

        int m_inputs, m_values;
        DetectorGeo* det_;
        std::vector<Track*> active_tracks;
        std::vector<Cluster*> unique_clusters;
        std::vector<int> active_dofs;
        ThreadPool* pool_;

        Eigen::Vector3d orig_pos, orig_rot;

        AlignmentFunctor(int inputs, int values, std::string det_name, std::vector<Track*> tracks,
                          std::vector<int> dofs, ThreadPool* pool)
            : m_inputs(inputs), m_values(values), active_tracks(tracks), active_dofs(dofs), pool_(pool) {

            det_ = &GeometryManager::getInstance().getDetector(det_name);
            orig_pos = det_->position;
            orig_rot = det_->orientation;

            std::unordered_set<Cluster*> cl_set;
            for (auto* t : active_tracks) {
                for (auto& cl : t->getClusters()) {
                    if (cl->detectorID() == det_->name) cl_set.insert(cl.get());
                }
            }
            unique_clusters.assign(cl_set.begin(), cl_set.end());
        }

        int inputs() const { return m_inputs; }
        int values() const { return m_values; }

        int operator()(const InputType &x, ValueType &fvec) const {
            Eigen::Vector3d dpos(0,0,0), drot(0,0,0);
            for (size_t i = 0; i < (size_t)m_inputs; ++i) {
                int dof = active_dofs[i];
                if (dof == 0) dpos.x() = x(i);
                if (dof == 1) dpos.y() = x(i);
                if (dof == 2) dpos.z() = x(i);
                if (dof == 3) drot.x() = x(i);
                if (dof == 4) drot.y() = x(i);
                if (dof == 5) drot.z() = x(i);
            }

            det_->position = orig_pos + dpos;
            det_->orientation = orig_rot + drot;
            det_->updateMatrices(); 

            for (auto* c : unique_clusters) {
                double lx = det_->getLocalX(c->column());
                double ly = det_->getLocalY(c->row());
                c->setGlobal(det_->R * Eigen::Vector3d(lx, ly, 0.0) + det_->position);
            }

            // Refitting every track touching this detector is the hottest loop in
            // the alignment procedure - it runs for every LM function/Jacobian
            // evaluation, every iteration, every detector. Chunked across worker
            // threads the same way Tracking4D/MaskCreator do.
            size_t n = active_tracks.size();
            size_t n_threads = pool_ ? pool_->getThreadCount() : 0;
            if (pool_ && n_threads > 1 && n > n_threads) {
                size_t chunk = (n + n_threads - 1) / n_threads;
                std::vector<std::future<void>> futures;
                for (size_t start = 0; start < n; start += chunk) {
                    size_t end = std::min(start + chunk, n);
                    futures.push_back(pool_->submit([this, start, end, &fvec]() {
                        for (size_t t = start; t < end; ++t) {
                            Track* track = active_tracks[t];
                            fvec(t) = track->fit() ? std::sqrt(track->chi2()) : 1e6;
                        }
                    }));
                }
                for (auto& f : futures) f.get();
            } else {
                for (size_t t = 0; t < n; ++t) {
                    Track* track = active_tracks[t];
                    fvec(t) = track->fit() ? std::sqrt(track->chi2()) : 1e6;
                }
            }
            return 0;
        }
    };


    AlignmentTrackChi2::AlignmentTrackChi2(Configuration cfg, Configuration g_cfg, ThreadPool* pool)
        : Module(std::move(cfg), std::move(g_cfg), pool) {
        iterations_ = config.get<int>("iterations", 5);
        max_chi2ndof_ = config.get<double>("max_track_chi2ndof", 30.0);
        prune_tracks_ = config.get<bool>("prune_tracks", false); 
        align_position_ = config.get<bool>("align_position", true);
        align_orientation_ = config.get<bool>("align_orientation", true);
        type_filter_ = config.get<std::string>("type", "");
        // Matches corry's global `number_of_tracks` (see FrameworkManager):
        // one framework-wide track target, not a module-local setting. 0
        // means unlimited.
        max_tracks_ = global_config.get<size_t>("number_of_tracks", 0);
        fixed_planes_ = config.getArray<std::string>("fixed_planes");
    }

    void AlignmentTrackChi2::initialize() {
        WR_LOG(STATUS, "Initializing AlignmentTrackChi2 (" + std::to_string(iterations_) + " iterations)");
        auto& geo = GeometryManager::getInstance();
        
        for (const auto& name : geo.getDetectorNames()) {
            if (name == geo.getReferenceName() || !geo.getDetector(name).isPixelDetector()) continue;
            std::string dir = getName() + "/" + name;
            PlotManager::getInstance().registerGraph(dir, "alignment_correction_displacementX");
            PlotManager::getInstance().registerGraph(dir, "alignment_correction_displacementY");
            PlotManager::getInstance().registerGraph(dir, "alignment_correction_rotationX");
            PlotManager::getInstance().registerGraph(dir, "alignment_correction_rotationY");
            PlotManager::getInstance().registerGraph(dir, "alignment_correction_rotationZ");
        }
    }

    void AlignmentTrackChi2::run(DataBatch& batch) {
        std::lock_guard<std::mutex> lock(track_mtx_);
        // No module-level track-count cap here, matching corry: it relies on
        // the framework-level `number_of_tracks` stopping check (see
        // FrameworkManager) rather than this module enforcing its own exact
        // boundary. That check only runs after a full event, so the total
        // can overshoot by one event's worth of tracks.
        for (auto& track : batch.tracks) {
            total_evaluated_tracks_++;
            if (!prune_tracks_ || track->getChi2ndof() <= max_chi2ndof_) {
                global_tracks_.push_back(track);
            }
        }
    }

    void AlignmentTrackChi2::finalize() {
        auto& geo = GeometryManager::getInstance();
        std::string ref = geo.getReferenceName();
        if (global_tracks_.empty()) return;

        size_t discarded = total_evaluated_tracks_ - global_tracks_.size();
        WR_LOG(INFO, "Discarded " + std::to_string(discarded) + " input tracks.");
        WR_LOG(STATUS, "Starting iterative alignment (Eigen Levenberg-Marquardt) with " + std::to_string(global_tracks_.size()) + " tracks.");

        std::vector<int> active_dofs;
        if (align_position_) { active_dofs.push_back(0); active_dofs.push_back(1); }
        if (align_orientation_) { active_dofs.push_back(3); active_dofs.push_back(4); active_dofs.push_back(5); }
        if (active_dofs.empty()) return;

        for (int iter = 0; iter < iterations_; ++iter) {
            WR_LOG(INFO, "--- Alignment iteration " + std::to_string(iter + 1) + " of " + std::to_string(iterations_) + " ---");

            for (const auto& name : geo.getDetectorNames()) {
                auto& det = geo.getDetector(name);
                if (name == ref || !det.isPixelDetector()) continue;
                if (!type_filter_.empty() && det.type != type_filter_) continue;
                if (std::find(fixed_planes_.begin(), fixed_planes_.end(), name) != fixed_planes_.end()) continue;

                std::vector<Track*> active_tracks;
                for (auto& t : global_tracks_) if (t->hasDetector(name)) active_tracks.push_back(t.get());
                if (active_tracks.empty()) continue;

                AlignmentFunctor functor(active_dofs.size(), active_tracks.size(), name, active_tracks, active_dofs, thread_pool);
                Eigen::NumericalDiff<AlignmentFunctor> numDiff(functor);
                Eigen::LevenbergMarquardt<Eigen::NumericalDiff<AlignmentFunctor>> lm(numDiff);

                lm.parameters.maxfev = 1000;
                lm.parameters.ftol = 1e-9; 
                lm.parameters.xtol = 1e-7;

                Eigen::VectorXd x = Eigen::VectorXd::Zero(active_dofs.size());
                lm.minimize(x);
                
                Eigen::VectorXd fvec(active_tracks.size());
                functor(x, fvec); 

                Eigen::Vector3d dpos(0,0,0), drot(0,0,0);
                for (size_t i = 0; i < active_dofs.size(); ++i) {
                    if (active_dofs[i] == 0) dpos.x() = x(i);
                    if (active_dofs[i] == 1) dpos.y() = x(i);
                    if (active_dofs[i] == 3) drot.x() = x(i);
                    if (active_dofs[i] == 4) drot.y() = x(i);
                    if (active_dofs[i] == 5) drot.z() = x(i);
                }

                std::ostringstream oss;
                oss << name << "/" << iter << " dT(" 
                    << std::fixed << std::setprecision(3) << dpos.x() * 1000.0 << "um," 
                    << dpos.y() * 1000.0 << "um,0mm) dR("
                    << drot.x() * 180.0/M_PI << "deg," << drot.y() * 180.0/M_PI << "deg," << drot.z() * 180.0/M_PI << "deg)";
                WR_LOG(DEBUG, oss.str());

                std::string dir = getName() + "/" + name;
                PlotManager::getInstance().addPoint(dir, "alignment_correction_displacementX", iter, dpos.x() * 1000.0);
                PlotManager::getInstance().addPoint(dir, "alignment_correction_displacementY", iter, dpos.y() * 1000.0);
                PlotManager::getInstance().addPoint(dir, "alignment_correction_rotationX", iter, drot.x() * 180.0/M_PI);
                PlotManager::getInstance().addPoint(dir, "alignment_correction_rotationY", iter, drot.y() * 180.0/M_PI);
                PlotManager::getInstance().addPoint(dir, "alignment_correction_rotationZ", iter, drot.z() * 180.0/M_PI);
            }
        }

        WR_LOG(STATUS, "--- Final Alignment Parameters ---");
        for (const auto& name : geo.getDetectorNames()) {
            auto& det = geo.getDetector(name);
            if (name == ref || !det.isPixelDetector()) continue;

            std::ostringstream oss_r;
            oss_r << std::fixed << std::setprecision(4) << det.orientation(0)*180.0/M_PI << "deg, " << det.orientation(1)*180.0/M_PI << "deg, " << det.orientation(2)*180.0/M_PI << "deg";

            WR_LOG(STATUS, name + " new alignment:");
            WR_LOG(STATUS, "  T(" + Units::display(det.position(0)) + "," + Units::display(det.position(1)) + "," + Units::display(det.position(2)) + ")");
            WR_LOG(STATUS, "  R(" + oss_r.str() + ")");
        }

        std::string out_geo = global_config.get<std::string>("detectors_file_updated", "aligned.geo");
        geo.saveGeometry(out_geo);
        WR_LOG(STATUS, "Alignment finished. Updated geometry saved to: " + out_geo);
    }

    REGISTER_MODULE(AlignmentTrackChi2)
}