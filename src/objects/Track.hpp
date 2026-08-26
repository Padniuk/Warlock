#ifndef WARLOCK_TRACK_HPP
#define WARLOCK_TRACK_HPP

#include "Cluster.hpp"
#include "core/detector/GeometryManager.hpp"
#include <Eigen/Dense>
#include <vector>
#include <memory>
#include <string>
#include <map>
#include <cmath>

namespace framework {
    class Track : public Object {
    public:
        Track() : Object(0) {}
        virtual ~Track() = default;

        void addCluster(std::shared_ptr<Cluster> c) { clusters_.push_back(c); }
        const std::vector<std::shared_ptr<Cluster>>& getClusters() const { return clusters_; }

        void addAssociatedCluster(std::shared_ptr<Cluster> cluster) { associated_clusters_.push_back(cluster); }
        const std::vector<std::shared_ptr<Cluster>>& getAssociatedClusters() const { return associated_clusters_; }

        void setClosestCluster(std::shared_ptr<Cluster> cluster) { closest_cluster_ = cluster; }
        std::shared_ptr<Cluster> getClosestCluster() const { return closest_cluster_; }

        void setParams(Eigen::Vector4d p) { params_ = std::move(p); }
        const Eigen::Vector4d& params() const { return params_; }

        void setChi2(double c) { chi2_ = c; }
        double chi2() const { return chi2_; }
        
        void setNdof(int n) { ndof_ = n; }
        int ndof() const { return ndof_; }

        double getChi2ndof() const { return ndof_ > 0 ? (chi2_ / static_cast<double>(ndof_)) : -1.0; }
        bool isFitted() const { return is_fitted_; }

        Eigen::Vector2d positionAt(double z) const { 
            return {params_(0) * z + params_(1), params_(2) * z + params_(3)}; 
        }

        Eigen::Vector3d interceptAt(const DetectorGeo& det) const {
            if (det.orientation.norm() < 1e-9) {
                return {params_(0) * det.position.z() + params_(1), params_(2) * det.position.z() + params_(3), det.position.z()};
            }
            Eigen::Vector3d planeN = det.R * Eigen::Vector3d(0, 0, 1);
            
            Eigen::Vector3d trk_state(params_(1), params_(3), 0.0);
            Eigen::Vector3d trk_dir(params_(0), params_(2), 1.0);
            
            double pathLength = -(trk_state - det.position).dot(planeN) / trk_dir.dot(planeN);
            return trk_state + pathLength * trk_dir;
        }

        bool hasDetector(const std::string& name) const {
            for(auto const& c : clusters_) if(c->detectorID() == name) return true;
            for(auto const& c : associated_clusters_) if(c->detectorID() == name) return true;
            return false;
        }

        // Matches corry's Track::isIsolated(): true unless some module has
        // explicitly marked this track as too close to another track at one
        // of the planes it was evaluated on via markAsNotIsolatedAtPlane().
        bool isIsolated() const {
            for (auto const& [name, iso] : is_isolated_at_) if (!iso) return false;
            return true;
        }
        void markAsNotIsolatedAtPlane(const std::string& det_name) { is_isolated_at_[det_name] = false; }

        virtual bool fit(const std::string& shift_det = "", const Eigen::Vector3d& shift = Eigen::Vector3d::Zero()) {
            size_t n = clusters_.size();
            if (n < 2) return false;

            // Exact 2-point line through both clusters (no fitting needed).
            if (n == 2 && shift_det.empty()) {
                double z0 = clusters_[0]->global().z();
                double z1 = clusters_[1]->global().z();
                if (std::abs(z1 - z0) < 1e-6) {
                    return false;
                }

                double dz = z1 - z0;
                double slope_x = (clusters_[1]->global().x() - clusters_[0]->global().x()) / dz;
                double intercept_x = clusters_[0]->global().x() - slope_x * z0;

                double slope_y = (clusters_[1]->global().y() - clusters_[0]->global().y()) / dz;
                double intercept_y = clusters_[0]->global().y() - slope_y * z0;

                params_ = {slope_x, intercept_x, slope_y, intercept_y};
                chi2_ = 0.0;
                ndof_ = 0;
                is_fitted_ = true;
                return true;
            }

            // Analytic least-squares fit.
            bool is_coupled = false;
            for (auto& c : clusters_) {
                if (c->detector()->orientation.norm() > 1e-9) is_coupled = true;
            }

            if (!is_coupled) {
                double S_x = 0, Sz_x = 0, Szz_x = 0, Sx_x = 0, Sxz_x = 0;
                double S_y = 0, Sz_y = 0, Szz_y = 0, Sy_y = 0, Syz_y = 0;

                double z_mean = 0;
                for (size_t i = 0; i < n; ++i) z_mean += clusters_[i]->global().z();
                z_mean /= n;

                for (size_t i = 0; i < n; ++i) {
                    auto* det = clusters_[i]->detector();
                    double x = clusters_[i]->global().x(); 
                    double y = clusters_[i]->global().y(); 
                    double z = clusters_[i]->global().z();
                    
                    if (det->name == shift_det) { x += shift.x(); y += shift.y(); z += shift.z(); }

                    double z_prime = z - z_mean; 

                    // Diagonal weights from spatial_res_x/y.
                    double wx = 1.0 / (det->spatial_res_x * det->spatial_res_x);
                    double wy = 1.0 / (det->spatial_res_y * det->spatial_res_y);

                    S_x += wx;         Sz_x += wx * z_prime;      Szz_x += wx * z_prime * z_prime;
                    Sx_x += wx * x;    Sxz_x += wx * x * z_prime;

                    S_y += wy;         Sz_y += wy * z_prime;      Szz_y += wy * z_prime * z_prime;
                    Sy_y += wy * y;    Syz_y += wy * y * z_prime;
                }

                double det_x = S_x * Szz_x - Sz_x * Sz_x;
                double det_y = S_y * Szz_y - Sz_y * Sz_y;

                if (std::abs(det_x) < 1e-12 || std::abs(det_y) < 1e-12) {
                    return false;
                }

                double c_x_prime = (Sx_x * Szz_x - Sxz_x * Sz_x) / det_x;
                double m_x = (Sxz_x * S_x - Sx_x * Sz_x) / det_x;
                double c_x = c_x_prime - m_x * z_mean; 

                double c_y_prime = (Sy_y * Szz_y - Syz_y * Sz_y) / det_y;
                double m_y = (Syz_y * S_y - Sy_y * Sz_y) / det_y;
                double c_y = c_y_prime - m_y * z_mean; 

                params_ = {m_x, c_x, m_y, c_y};
            } else {
                Eigen::Matrix4d mat = Eigen::Matrix4d::Zero();
                Eigen::Vector4d vec = Eigen::Vector4d::Zero();

                for (size_t i = 0; i < n; ++i) {
                    auto* det = clusters_[i]->detector();
                    double x = clusters_[i]->global().x(); 
                    double y = clusters_[i]->global().y(); 
                    double z = clusters_[i]->global().z();
                    if (det->name == shift_det) { x += shift.x(); y += shift.y(); z += shift.z(); }

                    // Full covariance weights, rotated into the global frame.
                    Eigen::Matrix3d V_loc = Eigen::Matrix3d::Zero();
                    V_loc(0,0) = det->spatial_res_x * det->spatial_res_x;
                    V_loc(1,1) = det->spatial_res_y * det->spatial_res_y;
                    V_loc(2,2) = 1e-9;
                    Eigen::Matrix3d V_glob = det->R * V_loc * det->R_inv;
                    Eigen::Matrix2d V;
                    V << V_glob(0,0), V_glob(0,1), V_glob(1,0), V_glob(1,1);
                    Eigen::Matrix2d V_i = V.inverse();

                    double v00 = V_i(0,0), v01 = V_i(0,1);
                    double v10 = V_i(1,0), v11 = V_i(1,1);

                    double vx = v00 * x + v01 * y;
                    double vy = v10 * x + v11 * y;
                    
                    vec(0) += vx; vec(1) += vx * z; vec(2) += vy; vec(3) += vy * z;

                    mat(0,0) += v00;       mat(0,1) += v00 * z;   mat(0,2) += v01;       mat(0,3) += v01 * z;
                    mat(1,0) += v00 * z;   mat(1,1) += v00 * z*z; mat(1,2) += v01 * z;   mat(1,3) += v01 * z*z;
                    mat(2,0) += v10;       mat(2,1) += v10 * z;   mat(2,2) += v11;       mat(2,3) += v11 * z;
                    mat(3,0) += v10 * z;   mat(3,1) += v10 * z*z; mat(3,2) += v11 * z;   mat(3,3) += v11 * z*z;
                }
                Eigen::Vector4d res = mat.colPivHouseholderQr().solve(vec);
                params_ = {res(1), res(0), res(3), res(2)};
            }

            if (std::isnan(params_(0)) || std::isnan(params_(1)) || std::isnan(params_(2)) || std::isnan(params_(3)) ||
                std::isinf(params_(0)) || std::isinf(params_(1)) || std::isinf(params_(2)) || std::isinf(params_(3))) {
                return false;
            }

            chi2_ = 0;
            for (size_t i = 0; i < n; ++i) {
                auto* det = clusters_[i]->detector();
                Eigen::Vector3d global_pred = interceptAt(*det);
                Eigen::Vector3d hit_pos = clusters_[i]->global();
                if (det->name == shift_det) hit_pos += shift;

                if (det->orientation.norm() < 1e-9) {
                    chi2_ += std::pow(hit_pos.x() - global_pred.x(), 2) / (det->spatial_res_x * det->spatial_res_x) +
                             std::pow(hit_pos.y() - global_pred.y(), 2) / (det->spatial_res_y * det->spatial_res_y);
                } else {
                    Eigen::Vector3d local_pred = det->R_inv * (global_pred - det->position);
                    Eigen::Vector3d local_hit = det->R_inv * (hit_pos - det->position);
                    chi2_ += std::pow(local_hit.x() - local_pred.x(), 2) / (det->spatial_res_x * det->spatial_res_x) +
                             std::pow(local_hit.y() - local_pred.y(), 2) / (det->spatial_res_y * det->spatial_res_y);
                }
            }
            
            if (std::isnan(chi2_) || std::isinf(chi2_)) {
                return false;
            }

            ndof_ = static_cast<int>(2 * n - 4);
            is_fitted_ = true;
            return true;
        }

    protected:
        std::vector<std::shared_ptr<Cluster>> clusters_;
        std::vector<std::shared_ptr<Cluster>> associated_clusters_;
        std::shared_ptr<Cluster> closest_cluster_;
        std::map<std::string, bool> is_isolated_at_;

        Eigen::Vector4d params_{0,0,0,0}; 
        double chi2_{0.0};
        int ndof_{0};
        bool is_fitted_{false};
    };
}
#endif