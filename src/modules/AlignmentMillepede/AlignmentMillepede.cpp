/** @file AlignmentMillepede.cpp */
#include "AlignmentMillepede.hpp"
#include "core/ModuleFactory.hpp"
#include "core/detector/GeometryManager.hpp"
#include "core/utils/PlotManager.hpp"
#include "core/utils/Logger.hpp"
#include "core/utils/Units.hpp"
#include <algorithm>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <future>
#include <atomic>

namespace framework {

    AlignmentMillepede::AlignmentMillepede(Configuration cfg, Configuration g_cfg, ThreadPool* pool)
        : Module(std::move(cfg), std::move(g_cfg), pool) {
        
        iterations_ = config.get<int>("iterations", 5);
        // Accepts either detector names or detector types, as a single
        // bare value or a list (see Configuration::getArray()) - matches
        // Tracking4D::exclude_detectors_, so a plane excluded/included from
        // the track fit and this alignment can be kept in sync from config
        // alone, without touching role in the geometry file.
        exclude_detectors_ = config.getArray<std::string>("exclude_detectors");
        dofs_ = config.getArray<bool>("dofs");
        if (dofs_.size() < 6) dofs_ = {true, true, false, true, true, true};
        sigmas_ = config.getArray<double>("sigmas");
        if (sigmas_.size() < 6) sigmas_ = {0.05, 0.05, 0.5, 0.005, 0.005, 0.005};
        
        // Matches corry's global `number_of_tracks` (see FrameworkManager):
        // one framework-wide track target, not a module-local setting. 0
        // means unlimited.
        max_tracks_ = global_config.get<size_t>("number_of_tracks", 0);
        residual_cut_ = config.get<double>("residual_cut", 0.05);
        residual_cut_init_ = config.get<double>("residual_cut_init", 0.6);
        number_of_stddev_ = config.get<int>("number_of_stddev", 0);
        convergence_cut_ = config.get<double>("convergence", 0.00001);
    }

    bool AlignmentMillepede::isExcluded(const DetectorGeo& det) const {
        return std::find(exclude_detectors_.begin(), exclude_detectors_.end(), det.name) != exclude_detectors_.end() ||
               std::find(exclude_detectors_.begin(), exclude_detectors_.end(), det.type) != exclude_detectors_.end();
    }

    void AlignmentMillepede::initialize() {
        WR_LOG(STATUS, "Initializing AlignmentMillepede...");
        auto& geo = GeometryManager::getInstance();

        unsigned int index = 0;
        for (const auto& name : geo.getDetectorNames()) {
            auto& det = geo.getDetector(name);
            if (!det.isPixelDetector() || det.role == "auxiliary") continue;
            if (isExcluded(det)) continue;
            m_millePlanes[name] = index;
            ++index;
            std::string dir = getName() + "/" + name;
            PlotManager::getInstance().registerGraph(dir, "alignment_correction_displacementX");
            PlotManager::getInstance().registerGraph(dir, "alignment_correction_displacementY");
            PlotManager::getInstance().registerGraph(dir, "alignment_correction_rotationX");
            PlotManager::getInstance().registerGraph(dir, "alignment_correction_rotationY");
            PlotManager::getInstance().registerGraph(dir, "alignment_correction_rotationZ");
        }

        WR_LOG(INFO, "Using the following degrees of freedom:");
        const std::vector<std::string> labels = {"Translation X", "Translation Y", "Translation Z", "Rotation X", "Rotation Y", "Rotation Z"};
        for(unsigned int i = 0; i < 6; ++i) WR_LOG(INFO, labels[i] + "\t" + (dofs_[i] ? "ON" : "OFF"));
    }

    void AlignmentMillepede::run(DataBatch& batch) {
        std::lock_guard<std::mutex> lock(track_mtx_);
        // No module-level track-count cap here, matching corry: it relies on
        // the framework-level `number_of_tracks` stopping check (see
        // FrameworkManager) rather than this module enforcing its own exact
        // boundary. That check only runs after a full event, so the total
        // can overshoot by one event's worth of tracks.
        for (auto& track : batch.tracks) {
            total_evaluated_tracks_++;
            global_tracks_.push_back(track);
        }
    }

    void AlignmentMillepede::finalize() {
        size_t nPlanes = m_millePlanes.size();
        WR_LOG(STATUS, "Millepede alignment");
        WR_LOG(INFO, "Aligning " + std::to_string(nPlanes) + " planes");

        const size_t nParameters = 6 * nPlanes;
        for(int iteration = 0; iteration < iterations_; ++iteration) {
            setConstraints(nPlanes);
            reset(nPlanes, 100.);
            
            WR_LOG(INFO, "Feeding Millepede with " + std::to_string(global_tracks_.size()) + " tracks...");
            
            m_equations.assign(global_tracks_.size(), std::vector<Equation>());
            std::atomic<unsigned int> nSkipped{0};
            std::atomic<unsigned int> nOutliers{0};
            std::mutex mtx;

            size_t n_threads = thread_pool->getThreadCount();
            size_t chunk = (global_tracks_.size() + n_threads - 1) / n_threads;
            std::vector<std::future<void>> futures;

            for (size_t i = 0; i < global_tracks_.size(); i += chunk) {
                size_t end = std::min(i + chunk, global_tracks_.size());
                futures.push_back(thread_pool->submit([&, i, end]() {
                    ThreadWorkspace ws(m_nagb, m_nalc);
                    unsigned int l_nSkipped = 0, l_nOutliers = 0;

                    for (size_t t = i; t < end; ++t) {
                        auto* track = global_tracks_[t].get();
                        size_t valid_hits = 0;
                        for (auto& c : track->getClusters()) if (m_millePlanes.count(c->detectorID())) valid_hits++;
                        
                        if(valid_hits != nPlanes) { l_nSkipped++; continue; }
                        
                        std::vector<Equation> local_eqs;
                        if(!putTrack(track, nPlanes, local_eqs, ws)) {
                            l_nOutliers++;
                        } else {
                            m_equations[t] = std::move(local_eqs);
                        }
                    }
                    nSkipped += l_nSkipped;
                    nOutliers += l_nOutliers;

                    std::lock_guard<std::mutex> lock(mtx);
                    for(size_t k = 0; k < m_nagb; ++k) {
                        m_bgvec[k] += ws.bgvec[k];
                        for(size_t l = 0; l < m_nagb; ++l) m_cgmat[k][l] += ws.cgmat[k][l];
                    }
                }));
            }
            for (auto& f : futures) f.get();

            if(nSkipped > 0) WR_LOG(INFO, "Skipped " + std::to_string(nSkipped) + " tracks with less than " + std::to_string(nPlanes) + " clusters.");
            if(nOutliers > 0) WR_LOG(INFO, "Rejected " + std::to_string(nOutliers) + " outlier tracks.");

            WR_LOG(INFO, "Determining global parameters...");
            if(!fitGlobal()) { WR_LOG(ERROR, "Global fit failed."); break; }

            double converg = 0.;
            for(size_t i = 0; i < nParameters; ++i) converg += std::abs(m_dparm[i]);
            converg /= static_cast<double>(nParameters);
            
            std::ostringstream oss;
            oss << "Convergence: " << converg;
            WR_LOG(INFO, oss.str());
            WR_LOG(INFO, "Updating geometry...");
            updateGeometry();

            futures.clear();
            for (size_t i = 0; i < global_tracks_.size(); i += chunk) {
                size_t end = std::min(i + chunk, global_tracks_.size());
                futures.push_back(thread_pool->submit([&, i, end]() {
                    auto& geo = GeometryManager::getInstance();
                    for (size_t t = i; t < end; ++t) {
                        auto* track = global_tracks_[t].get();
                        for (auto& c : track->getClusters()) {
                            if (m_millePlanes.count(c->detectorID())) {
                                auto& det = geo.getDetector(c->detectorID());
                                // det.R was just refreshed by updateGeometry()'s
                                // updateMatrices() call, so it can be reused
                                // here instead of rebuilding it per cluster.
                                double lx = det.getLocalX(c->column());
                                double ly = det.getLocalY(c->row());
                                c->setGlobal(det.R * Eigen::Vector3d(lx, ly, 0.0) + det.position);
                            }
                        }
                        track->fit();
                    }
                }));
            }
            for (auto& f : futures) f.get();

            for(const auto& [name, plane] : m_millePlanes) {
                std::string dir = getName() + "/" + name;
                PlotManager::getInstance().addPoint(dir, "alignment_correction_displacementX", iteration, m_dparm[plane + 0 * nPlanes] * 1000.0);
                PlotManager::getInstance().addPoint(dir, "alignment_correction_displacementY", iteration, m_dparm[plane + 1 * nPlanes] * 1000.0);
                PlotManager::getInstance().addPoint(dir, "alignment_correction_rotationX", iteration, m_dparm[plane + 3 * nPlanes] * (180.0/M_PI));
                PlotManager::getInstance().addPoint(dir, "alignment_correction_rotationY", iteration, m_dparm[plane + 4 * nPlanes] * (180.0/M_PI));
                PlotManager::getInstance().addPoint(dir, "alignment_correction_rotationZ", iteration, m_dparm[plane + 5 * nPlanes] * (180.0/M_PI));
            }

            if(converg < convergence_cut_) break;
        }

        auto& geo = GeometryManager::getInstance();
        std::string out_geo = global_config.get<std::string>("detectors_file_updated", "mille_aligned.geo");
        geo.saveGeometry(out_geo);
        WR_LOG(STATUS, "Millepede Alignment Finished. Geometry Saved.");
    }

    void AlignmentMillepede::setConstraints(const size_t nPlanes) {
        auto& geo = GeometryManager::getInstance();
        double avgz = 0., varz = 0.0;
        for(const auto& [name, idx] : m_millePlanes) avgz += geo.getDetector(name).position.z();
        avgz /= static_cast<double>(nPlanes);
        for(const auto& [name, idx] : m_millePlanes) {
            double dz = geo.getDetector(name).position.z() - avgz; varz += dz * dz;
        }
        varz /= static_cast<double>(nPlanes);

        const size_t nParameters = 6 * nPlanes;
        std::vector<double> ftx(nParameters, 0.), fty(nParameters, 0.), ftz(nParameters, 0.);
        std::vector<double> frx(nParameters, 0.), fry(nParameters, 0.), frz(nParameters, 0.);
        std::vector<double> fscaz(nParameters, 0.), shearx(nParameters, 0.), sheary(nParameters, 0.);

        m_constraints.clear();
        for(const auto& [name, i] : m_millePlanes) {
            const double sz = (geo.getDetector(name).position.z() - avgz) / varz;
            ftx[i] = 1.0; fty[i + nPlanes] = 1.0; ftz[i + 2 * nPlanes] = 1.0;
            frx[i + 3 * nPlanes] = i >= 4 ? 1.0 : -1.0;
            fry[i + 4 * nPlanes] = i >= 4 ? 1.0 : -1.0;
            frz[i + 5 * nPlanes] = 1.0;
            shearx[i] = sz; sheary[i + nPlanes] = sz; fscaz[i + 2 * nPlanes] = sz;
        }

        // Matches corry's own hardcoded constraint-activation table: the
        // Z-translation constraint (ftz) and the Z-shear/scale constraint
        // (fscaz) are permanently disabled regardless of the dofs config -
        // only ftx/shearx/fty/sheary/frx/fry/frz are ever added. Z-translation
        // and the tilt angles (rotation X/Y) are a nearly degenerate direction
        // in track-based alignment, so gating ftz on dofs_[2] would change how
        // the fit apportions correction between shifting a plane in Z and
        // tilting it, diverging rotation X/Y from corry's output.
        if(dofs_[0]) addConstraint(ftx, 0.0);
        if(dofs_[0]) addConstraint(shearx, 0.);
        if(dofs_[1]) addConstraint(fty, 0.0);
        if(dofs_[1]) addConstraint(sheary, 0.);
        if(dofs_[3]) addConstraint(frx, 0.0);
        if(dofs_[4]) addConstraint(fry, 0.0);
        if(dofs_[5]) addConstraint(frz, 0.0);
    }

    void AlignmentMillepede::addConstraint(const std::vector<double>& dercs, const double rhs) {
        Constraint c; c.rhs = rhs; c.coefficients = dercs;
        m_constraints.push_back(c);
    }

    bool AlignmentMillepede::putTrack(Track* track, const size_t nPlanes, std::vector<Equation>& equations, ThreadWorkspace& ws) {
        auto& geo = GeometryManager::getInstance();
        track->fit();

        // Track intercept at the first plane's own Z, matching corry's
        // exact 3D intercept convention (not a straight-line extrapolation
        // to Z=0).
        auto& det0 = geo.getDetector(track->getClusters().front()->detectorID());
        Eigen::Vector3d intercept0 = track->interceptAt(det0);
        const double tx = intercept0.x(); 
        const double ty = intercept0.y(); 

        for(auto& cluster : track->getClusters()) {
            auto plane_it = m_millePlanes.find(cluster->detectorID());
            if(plane_it == m_millePlanes.end()) continue;
            auto& detector = geo.getDetector(cluster->detectorID());

            // detector.R is kept in sync with detector.orientation by
            // updateGeometry() (and by the initial geometry load), so it can
            // be used directly instead of rebuilding the same rotation matrix
            // from scratch for every cluster of every track.
            const Eigen::Matrix3d& R = detector.R;
            double nx = R(0,2) / R(2,2);
            double ny = R(1,2) / R(2,2);

            const double xg = cluster->global().x();
            const double yg = cluster->global().y();
            const double zg = cluster->global().z();
            const double zl = zg - detector.position.z();
            const double xl = xg - detector.position.x();
            const double yl = yg - detector.position.y();

            std::vector<double> nonlinear = {nx, ny, 1., -yl + zl * ny, -xl + zl * nx, xl * ny - yl * nx};
            const double den = 1. + tx * nx + ty * ny;
            for(auto& a : nonlinear) a /= den;

            const double errx = detector.spatial_res_x;
            const double erry = detector.spatial_res_y;
            const unsigned int plane = plane_it->second;

            // ws.dergb/dernl/dernli/dernls need no explicit reset here:
            // addEquation() self-clears every dergb entry it reads (a full scan
            // it already performs), and dernl/dernli/dernls are only ever read
            // at indices where dergb is non-zero - which are always freshly
            // overwritten below before that read happens.
            std::vector<double> derlc = {1., zg, 0., 0.};

            if(dofs_[0]) ws.dergb[plane] = -1.;
            if(dofs_[4]) ws.dergb[4 * nPlanes + plane] = -zl;
            if(dofs_[5]) ws.dergb[5 * nPlanes + plane] = yl;
            for(unsigned int i = 0; i < 6; ++i) {
                if(!dofs_[i]) continue;
                ws.dergb[i * nPlanes + plane] += tx * nonlinear[i];
                ws.dernl[i * nPlanes + plane] = nonlinear[i];
                ws.dernli[i * nPlanes + plane] = 1;
                ws.dernls[i * nPlanes + plane] = tx;
            }
            addEquation(equations, derlc, ws, xg, errx);

            derlc = {0., 0., 1., zg};

            if(dofs_[1]) ws.dergb[nPlanes + plane] = -1.;
            if(dofs_[3]) ws.dergb[3 * nPlanes + plane] = -zl;
            if(dofs_[5]) ws.dergb[5 * nPlanes + plane] = -xl;
            for(unsigned int i = 0; i < 6; ++i) {
                if(!dofs_[i]) continue;
                ws.dergb[i * nPlanes + plane] += ty * nonlinear[i];
                ws.dernl[i * nPlanes + plane] = nonlinear[i];
                ws.dernli[i * nPlanes + plane] = 3;
                ws.dernls[i * nPlanes + plane] = ty;
            }
            addEquation(equations, derlc, ws, yg, erry);
        }

        std::vector<double> trackParams(2 * m_nalc + 2, 0.);
        return fitTrackLocal(equations, trackParams, false, 1, ws);
    }

    void AlignmentMillepede::addEquation(std::vector<Equation>& equations, const std::vector<double>& derlc, ThreadWorkspace& ws, const double rmeas, const double sigma) {
        if(sigma <= 0.) {
            // Caller relies on this scan to reset ws.dergb for the next
            // equation; if we bail out before it runs, fall back to an
            // explicit clear so the workspace invariant still holds.
            std::fill(ws.dergb.begin(), ws.dergb.end(), 0.);
            return;
        }
        Equation equation;
        equation.rmeas = rmeas;
        equation.weight = 1. / (sigma * sigma);
        for(unsigned int i = 0; i < m_nalc; ++i) {
            if(derlc[i] == 0.) continue;
            equation.indL.push_back(i);
            equation.derL.push_back(derlc[i]);
        }
        // This scan doubles as the reset of ws.dergb for the next call - every
        // entry is visited, and any non-zero one is zeroed right after being
        // read, so no separate std::fill over ws.dergb is needed by callers.
        for(size_t i = 0; i < ws.dergb.size(); ++i) {
            if(ws.dergb[i] == 0.) continue;
            equation.indG.push_back(i);
            equation.derG.push_back(ws.dergb[i]);
            equation.derNL.push_back(ws.dernl[i]);
            equation.indNL.push_back(ws.dernli[i]);
            equation.slopes.push_back(ws.dernls[i]);
            ws.dergb[i] = 0.;
        }
        equations.push_back(std::move(equation));
    }

    bool AlignmentMillepede::fitTrackLocal(const std::vector<Equation>& equations, std::vector<double>& trackParams, const bool singlefit, const unsigned int iteration, ThreadWorkspace& ws) {
        std::vector<double> blvec(m_nalc, 0.);
        std::vector<std::vector<double>> clmat(m_nalc, std::vector<double>(m_nalc, 0.));

        for(const auto& equation : equations) {
            double rmeas = equation.rmeas;
            for(size_t i = 0; i < equation.derG.size(); ++i) rmeas -= equation.derG[i] * m_dparm[equation.indG[i]];
            const double w = equation.weight;
            for(size_t i = 0; i < equation.derL.size(); ++i) {
                const size_t j = equation.indL[i];
                blvec[j] += w * rmeas * equation.derL[i];
                for(unsigned int k = 0; k <= i; ++k) clmat[j][equation.indL[k]] += w * equation.derL[i] * equation.derL[k];
            }
        }

        int rank = invertMatrixLocal(clmat, blvec, m_nalc);
        for(unsigned int i = 0; i < m_nalc; ++i) {
            trackParams[2 * i] = blvec[i];
            trackParams[2 * i + 1] = std::sqrt(std::abs(clmat[i][i]));
        }

        double chi2 = 0.0;
        int ndf = 0;
        for(const auto& equation : equations) {
            double rmeas = equation.rmeas;
            for(size_t i = 0; i < equation.derL.size(); ++i) rmeas -= equation.derL[i] * blvec[equation.indL[i]];
            for(size_t i = 0; i < equation.derG.size(); ++i) rmeas -= equation.derG[i] * m_dparm[equation.indG[i]];
            
            const double rescut = iteration <= 1 ? residual_cut_init_ : residual_cut_;
            if(std::abs(rmeas) >= rescut) return false;
            chi2 += equation.weight * rmeas * rmeas;
            ++ndf;
        }
        ndf -= rank;
        if(singlefit) return true;

        if(number_of_stddev_ != 0 && ndf > 0) {
            const double cut = chi2Limit(number_of_stddev_, ndf) * m_cfactr;
            if(chi2 > cut) return false;
        }

        trackParams[2 * m_nalc] = ndf;
        trackParams[2 * m_nalc + 1] = chi2;

        for(auto& r : ws.clcmat) std::fill(r.begin(), r.end(), 0.);
        unsigned int nagbn = 0;
        std::vector<int> indnz(m_nagb, -1);
        std::vector<int> indbk(m_nagb, 0);

        for(const auto& equation : equations) {
            double rmeas = equation.rmeas;
            for(size_t i = 0; i < equation.derG.size(); ++i) rmeas -= equation.derG[i] * m_dparm[equation.indG[i]];
            const double w = equation.weight;
            
            for(size_t i = 0; i < equation.derG.size(); ++i) {
                const size_t j = equation.indG[i];
                ws.bgvec[j] += w * rmeas * equation.derG[i];
                for(size_t k = 0; k < equation.derG.size(); ++k) ws.cgmat[j][equation.indG[k]] += w * equation.derG[i] * equation.derG[k];
            }

            for(size_t i = 0; i < equation.derG.size(); ++i) {
                const size_t j = equation.indG[i];
                int ik = indnz[j];
                if(ik == -1) {
                    indnz[j] = nagbn; indbk[nagbn] = j;
                    ik = nagbn; ++nagbn;
                }
                for(size_t k = 0; k < equation.derL.size(); ++k) ws.clcmat[ik][equation.indL[k]] += w * equation.derG[i] * equation.derL[k];
            }
        }

        multiplyAVAt(clmat, ws.clcmat, ws.corrm, m_nalc, nagbn);
        multiplyAX(ws.clcmat, blvec, ws.corrv, m_nalc, nagbn);
        for(size_t i = 0; i < nagbn; ++i) {
            const size_t j = indbk[i];
            ws.bgvec[j] -= ws.corrv[i];
            for(size_t k = 0; k < nagbn; ++k) ws.cgmat[j][indbk[k]] -= ws.corrm[i][k];
        }
        return true;
    }

    void AlignmentMillepede::updateGeometry() {
        auto& geo = GeometryManager::getInstance();
        size_t nPlanes = m_millePlanes.size();

        for(const auto& [name, plane] : m_millePlanes) {
            auto& det = geo.getDetector(name);
            det.position(0) += m_dparm[plane + 0 * nPlanes];
            det.position(1) += m_dparm[plane + 1 * nPlanes];
            det.position(2) += m_dparm[plane + 2 * nPlanes];
            det.orientation(0) += m_dparm[plane + 3 * nPlanes];
            det.orientation(1) += m_dparm[plane + 4 * nPlanes];
            det.orientation(2) += m_dparm[plane + 5 * nPlanes];
            // Keep det.R/R_inv/V_inv in sync with the just-updated position/
            // orientation. Without this, they stay frozen at whatever they
            // were before this alignment ran, so anything reading det.R
            // directly (e.g. Track::fit()'s rotated-plane path) would use a
            // stale rotation while cluster global positions already reflect
            // the new one.
            det.updateMatrices();
        }
    }

    bool AlignmentMillepede::reset(const size_t nPlanes, const double startfact) {
        m_nagb = 6 * nPlanes;
        const unsigned int nRows = m_nagb + m_constraints.size();
        
        m_bgvec.assign(nRows, 0.);
        m_cgmat.assign(nRows, std::vector<double>(nRows, 0.));
        m_dparm.assign(m_nagb, 0.);
        m_fixed.assign(m_nagb, true);
        m_psigm.assign(m_nagb, 0.);

        for(unsigned int i = 0; i < 6; ++i) {
            if(!dofs_[i]) continue;
            for(size_t j = i * nPlanes; j < (i + 1) * nPlanes; ++j) {
                m_fixed[j] = false; m_psigm[j] = sigmas_[i];
            }
        }

        m_cfactref = 1.; m_cfactr = std::max(1., startfact);
        return true;
    }

    bool AlignmentMillepede::fitGlobal() {
        m_diag.assign(m_nagb, 0.);
        std::vector<double> bgvecPrev(m_nagb, 0.);
        const size_t nTracks = m_equations.size();
        std::vector<std::vector<double>> localParams(nTracks, std::vector<double>(m_nalc, 0.));

        const unsigned int nMaxIterations = 10;
        unsigned int iteration = 1;
        unsigned int nGoodTracks = nTracks;

        while(iteration <= nMaxIterations) {
            if(nGoodTracks == 0) return false;
            WR_LOG(DEBUG, "Iteration " + std::to_string(iteration) + " (using " + std::to_string(nGoodTracks) + " tracks)");

            for(unsigned int i = 0; i < m_nagb; ++i) m_diag[i] = m_cgmat[i][i];

            unsigned int nFixed = 0;
            for(unsigned int i = 0; i < m_nagb; ++i) {
                if(m_fixed[i]) {
                    ++nFixed;
                    for(unsigned int j = 0; j < m_nagb; ++j) { m_cgmat[i][j] = 0.0; m_cgmat[j][i] = 0.0; }
                } else m_cgmat[i][i] += 1. / (m_psigm[i] * m_psigm[i]);
            }

            size_t nRows = m_nagb;
            for(const auto& constraint : m_constraints) {
                double sum = constraint.rhs;
                for(unsigned int j = 0; j < m_nagb; ++j) {
                    if(m_psigm[j] == 0.) {
                        m_cgmat[nRows][j] = 0.0; m_cgmat[j][nRows] = 0.0;
                    } else {
                        m_cgmat[nRows][j] = nGoodTracks * constraint.coefficients[j];
                        m_cgmat[j][nRows] = m_cgmat[nRows][j];
                    }
                    sum -= constraint.coefficients[j] * m_dparm[j];
                }
                m_cgmat[nRows][nRows] = 0.0;
                m_bgvec[nRows] = nGoodTracks * sum;
                ++nRows;
            }

            const int rank = invertMatrix(m_cgmat, m_bgvec, nRows);
            for(unsigned int i = 0; i < m_nagb; ++i) {
                m_dparm[i] += m_bgvec[i];
                bgvecPrev[i] = m_bgvec[i];
            }

            if(nRows - nFixed - static_cast<unsigned int>(rank) != 0) {
                WR_LOG(WARNING, "The rank defect of the symmetric " + std::to_string(nRows) + " by " +
                                     std::to_string(nRows) + " matrix is " +
                                     std::to_string(nRows - nFixed - static_cast<unsigned int>(rank)));
            }

            ++iteration;
            if(iteration == nMaxIterations) break;

            const double newcfactr = std::sqrt(m_cfactr);
            m_cfactr = (newcfactr > 1.2 * m_cfactref) ? newcfactr : m_cfactref;

            for(unsigned int i = 0; i < nRows; ++i) {
                m_bgvec[i] = 0.0;
                for(unsigned int j = 0; j < nRows; ++j) m_cgmat[i][j] = 0.0;
            }

            double chi2 = 0.; double ndof = 0.; nGoodTracks = 0;
            std::mutex mtx;
            size_t n_threads = thread_pool->getThreadCount();
            size_t chunk = (nTracks + n_threads - 1) / n_threads;
            std::vector<std::future<void>> futures;

            for (size_t i = 0; i < nTracks; i += chunk) {
                size_t end = std::min(i + chunk, nTracks);
                futures.push_back(thread_pool->submit([&, i, end, iteration]() {
                    ThreadWorkspace ws(m_nagb, m_nalc);
                    double l_chi2 = 0, l_ndof = 0; unsigned int l_nGoodTracks = 0;
                    std::vector<double> l_trackParams(2 * m_nalc + 2, 0.);

                    for (size_t t = i; t < end; ++t) {
                        if(m_equations[t].empty()) continue;
                        // Nonlinear derivative corrections are intentionally not
                        // reapplied each inner iteration: corry's own inner-refit
                        // loop has an equivalent no-op (a by-value loop variable
                        // silently discards the update), so its output always
                        // reuses the initial linearization from putTrack(). Matching
                        // that here keeps output aligned with corry - applying the
                        // correction measurably shifts the (weakly-constrained)
                        // rotation parameters relative to corry's output.
                        auto& equations = m_equations[t];
                        std::fill(l_trackParams.begin(), l_trackParams.end(), 0.);
                        bool ok = fitTrackLocal(equations, l_trackParams, false, iteration, ws);
                        for(unsigned int j = 0; j < m_nalc; ++j) localParams[t][j] = l_trackParams[2 * j];

                        if(ok) {
                            l_chi2 += l_trackParams[2 * m_nalc + 1];
                            l_ndof += l_trackParams[2 * m_nalc];
                            ++l_nGoodTracks;
                        } else {
                            m_equations[t].clear();
                        }
                    }
                    std::lock_guard<std::mutex> lock(mtx);
                    chi2 += l_chi2; ndof += l_ndof; nGoodTracks += l_nGoodTracks;
                    for(size_t k = 0; k < m_nagb; ++k) {
                        m_bgvec[k] += ws.bgvec[k];
                        for(size_t l = 0; l < m_nagb; ++l) m_cgmat[k][l] += ws.cgmat[k][l];
                    }
                }));
            }
            for (auto& f : futures) f.get();

            std::ostringstream oss;
            oss << "Chi2 / DOF after re-fit: " << (chi2 / (ndof - static_cast<double>(nRows)));
            WR_LOG(DEBUG, oss.str());
        }

        printResults();
        return true;
    }

    int AlignmentMillepede::invertMatrix(std::vector<std::vector<double>>& v, std::vector<double>& b, const size_t n) {
        int rank = 0;
        const double eps = 0.0000000000001;
        std::vector<bool> flag(n, true), used_param(n, true);
        std::vector<double> r(n, 0.), c(n, 0.);

        for(size_t i = 0; i < n; i++) for(size_t j = 0; j <= i; j++) v[j][i] = v[i][j];
        for(size_t i = 0; i < n; i++) {
            for(size_t j = 0; j < n; j++) {
                if(std::abs(v[i][j]) >= r[i]) r[i] = std::abs(v[i][j]);
                if(std::abs(v[j][i]) >= c[i]) c[i] = std::abs(v[j][i]);
            }
        }
        for(size_t i = 0; i < n; i++) {
            if(r[i] != 0.0) r[i] = 1. / r[i];
            if(c[i] != 0.0) c[i] = 1. / c[i];
            if(eps >= r[i]) r[i] = 0.0;
            if(eps >= c[i]) c[i] = 0.0;
        }
        for(size_t i = 0; i < n; i++) for(size_t j = 0; j < n; j++) v[i][j] = std::sqrt(r[i]) * v[i][j] * std::sqrt(c[j]);
        for(size_t i = 0; i < n; i++) if(r[i] == 0. && c[i] == 0.) { flag[i] = false; used_param[i] = false; }

        for(size_t i = 0; i < n; i++) {
            double vkk = 0.0; size_t k = 0; bool found_k = false;
            for(size_t j = 0; j < n; j++) {
                if(flag[j] && (std::abs(v[j][j]) > std::max(std::abs(vkk), eps))) { vkk = v[j][j]; k = j; found_k = true; }
            }
            if(found_k) {
                rank++; flag[k] = false; vkk = 1.0 / vkk; v[k][k] = -vkk;
                for(size_t j = 0; j < n; j++) for(size_t jj = 0; jj < n; jj++) {
                    if(j != k && jj != k && used_param[j] && used_param[jj]) v[j][jj] = v[j][jj] - vkk * v[j][k] * v[k][jj];
                }
                for(size_t j = 0; j < n; j++) {
                    if(j != k && used_param[j]) { v[j][k] *= vkk; v[k][j] *= vkk; }
                }
            } else {
                for(size_t j = 0; j < n; j++) if(flag[j]) {
                    b[j] = 0.0; for(size_t l = 0; l < n; l++) { v[j][l] = 0.0; v[l][j] = 0.0; }
                }
                break;
            }
        }
        for(size_t i = 0; i < n; i++) for(size_t j = 0; j < n; j++) v[i][j] = std::sqrt(c[i]) * v[i][j] * std::sqrt(r[j]);
        std::vector<double> temp(n, 0.);
        for(size_t j = 0; j < n; j++) for(size_t jj = 0; jj < n; jj++) { v[j][jj] = -v[j][jj]; temp[j] += v[j][jj] * b[jj]; }
        for(size_t j = 0; j < n; j++) b[j] = temp[j];
        return rank;
    }

    int AlignmentMillepede::invertMatrixLocal(std::vector<std::vector<double>>& v, std::vector<double>& b, const size_t n) {
        int rank = 0;
        const double eps = 0.0000000000001;
        std::vector<bool> flag(n, true);
        std::vector<double> diag(n, 0.);
        for(size_t i = 0; i < n; i++) {
            diag[i] = std::abs(v[i][i]);
            for(size_t j = 0; j <= i; j++) v[j][i] = v[i][j];
        }

        for(size_t i = 0; i < n; i++) {
            double vkk = 0.0; size_t k = 0; bool found_k = false;
            for(size_t j = 0; j < n; j++) {
                if(flag[j] && (std::abs(v[j][j]) > std::max(std::abs(vkk), eps * diag[j]))) { vkk = v[j][j]; k = j; found_k = true; }
            }
            if(found_k) {
                rank++; flag[k] = false; vkk = 1.0 / vkk; v[k][k] = -vkk;
                for(size_t j = 0; j < n; j++) for(size_t jj = 0; jj < n; jj++) {
                    if(j != k && jj != k) v[j][jj] = v[j][jj] - vkk * v[j][k] * v[k][jj];
                }
                for(size_t j = 0; j < n; j++) if(j != k) { v[j][k] *= vkk; v[k][j] *= vkk; }
            } else {
                for(size_t j = 0; j < n; j++) if(flag[j]) { b[j] = 0.0; for(size_t l=0; l<n; l++) v[j][l] = 0.0; }
                break;
            }
        }
        std::vector<double> temp(n, 0.);
        for(size_t j = 0; j < n; j++) for(size_t jj = 0; jj < n; jj++) { v[j][jj] = -v[j][jj]; temp[j] += v[j][jj] * b[jj]; }
        for(size_t j = 0; j < n; j++) b[j] = temp[j];
        return rank;
    }

    double AlignmentMillepede::chi2Limit(const int n, const int nd) const {
        constexpr double sn[3] = {0.47523, 1.690140, 2.782170};
        constexpr double table[3][30] = {
            {1.0000, 1.1479, 1.1753, 1.1798, 1.1775, 1.1730, 1.1680, 1.1630, 1.1581, 1.1536, 1.1493, 1.1454, 1.1417, 1.1383, 1.1351, 1.1321, 1.1293, 1.1266, 1.1242, 1.1218, 1.1196, 1.1175, 1.1155, 1.1136, 1.1119, 1.1101, 1.1085, 1.1070, 1.1055, 1.1040},
            {4.0000, 3.0900, 2.6750, 2.4290, 2.2628, 2.1415, 2.0481, 1.9736, 1.9124, 1.8610, 1.8171, 1.7791, 1.7457, 1.7161, 1.6897, 1.6658, 1.6442, 1.6246, 1.6065, 1.5899, 1.5745, 1.5603, 1.5470, 1.5346, 1.5230, 1.5120, 1.5017, 1.4920, 1.4829, 1.4742},
            {9.0000, 5.9146, 4.7184, 4.0628, 3.6410, 3.3436, 3.1209, 2.9468, 2.8063, 2.6902, 2.5922, 2.5082, 2.4352, 2.3711, 2.3143, 2.2635, 2.2178, 2.1764, 2.1386, 2.1040, 2.0722, 2.0428, 2.0155, 1.9901, 1.9665, 1.9443, 1.9235, 1.9040, 1.8855, 1.8681}
        };
        if(nd < 1) return 0.;
        const int m = std::max(1, std::min(n, 3));
        if(nd <= 30) return table[m - 1][nd - 1];
        return ((sn[m - 1] + std::sqrt((double)(2 * nd - 3))) * (sn[m - 1] + std::sqrt((double)(2 * nd - 3)))) / (double)(2 * nd - 2);
    }

    bool AlignmentMillepede::multiplyAX(const std::vector<std::vector<double>>& a, const std::vector<double>& x, std::vector<double>& y, const unsigned int n, const unsigned int m) {
        for(unsigned int i = 0; i < m; ++i) {
            y[i] = 0.0;
            for(unsigned int j = 0; j < n; ++j) y[i] += a[i][j] * x[j];
        }
        return true;
    }

    bool AlignmentMillepede::multiplyAVAt(const std::vector<std::vector<double>>& v, const std::vector<std::vector<double>>& a, std::vector<std::vector<double>>& w, const unsigned int n, const unsigned int m) {
        for(unsigned int i = 0; i < m; ++i) {
            for(unsigned int j = 0; j <= i; ++j) {
                w[i][j] = 0.0;
                for(unsigned int k = 0; k < n; ++k) for(unsigned int l = 0; l < n; ++l) w[i][j] += a[i][k] * v[k][l] * a[j][l];
                w[j][i] = w[i][j];
            }
        }
        return true;
    }

    void AlignmentMillepede::printResults() {
        const std::string line(65, '-');
        WR_LOG(DEBUG, line);
        WR_LOG(DEBUG, " Result of fit for global parameters");
        WR_LOG(DEBUG, line);
        WR_LOG(DEBUG, "   I  Difference    Last step      Error        Pull Global corr.");
        WR_LOG(DEBUG, line);

        for(unsigned int i = 0; i < m_nagb; ++i) {
            double err = std::sqrt(std::abs(m_cgmat[i][i]));
            if(m_cgmat[i][i] < 0.0) err = -err;

            if(std::abs(m_cgmat[i][i] * m_diag[i]) > 0) {
                const double pull = m_dparm[i] / std::sqrt(m_psigm[i] * m_psigm[i] - m_cgmat[i][i]);
                const double gcor = std::sqrt(std::abs(1.0 - 1.0 / (m_cgmat[i][i] * m_diag[i])));
                std::ostringstream oss;
                oss << std::setprecision(3) << std::scientific << std::setw(3) << i << "   " << std::setw(10) << m_dparm[i] << "   " << std::setw(10) << m_bgvec[i] << "   " << std::setw(8) << std::setprecision(2) << err << "   " << std::setw(9) << std::setprecision(2) << pull << "   " << std::setw(9) << gcor;
                WR_LOG(DEBUG, oss.str());
            } else {
                std::ostringstream oss;
                oss << std::setw(3) << i << "   " << std::setw(10) << "OFF" << "   " << std::setw(10) << "OFF" << "   " << std::setw(8) << "OFF" << "   " << std::setw(9) << "OFF" << "   " << std::setw(9) << "OFF";
                WR_LOG(DEBUG, oss.str());
            }
            if((i + 1) % (m_nagb / 6) == 0) WR_LOG(DEBUG, line);
        }
    }

    REGISTER_MODULE(AlignmentMillepede)
}