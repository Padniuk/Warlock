#include "MaskCreator.hpp"
#include "core/ModuleFactory.hpp"
#include "core/detector/GeometryManager.hpp"
#include "core/utils/PlotManager.hpp"
#include "core/utils/Logger.hpp"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <future>

namespace framework {

    MaskCreator::MaskCreator(Configuration cfg, Configuration g_cfg, ThreadPool* pool)
        : Module(std::move(cfg), std::move(g_cfg), pool) {

        type_filter_      = config.get<std::string>("type", "");
        method_           = config.get<std::string>("method", "frequency");
        frequency_cut_    = config.get<double>("frequency_cut", 50.0);

        // Local Density parameters (matches Corry's "Ignore" module logic)
        density_bandwidth_ = config.get<int>("density_bandwidth", 10);
        sigma_max_        = config.get<double>("sigma_above_avg_max", 5.0);
        rate_max_         = config.get<double>("rate_max", 1.0);

        mask_dead_pixels_ = config.get<bool>("mask_dead_pixels", false);
    }

    void MaskCreator::initialize() {
        auto& geo = GeometryManager::getInstance();
        WR_LOG(STATUS, "Initializing MaskCreator. Method: " + method_ + " Filter: [" + type_filter_ + "]");

        for (const auto& name : geo.getDetectorNames()) {
            auto& det = geo.getDetector(name);

            if (!type_filter_.empty() && det.type != type_filter_) continue;
            if (!det.isPixelDetector()) continue;

            occupancy_matrices_[name] = Eigen::MatrixXd::Zero(det.n_pixels_x, det.n_pixels_y);
            detector_mutexes_[name] = std::make_unique<std::mutex>();

            // Pre-fill with any pixels already masked via an existing mask_file,
            // matching corry's maskmap being pre-seeded from previously known
            // masks rather than this run silently dropping them.
            // GeometryManager::loadGeometry() already parsed mask_file into
            // det.masked_pixels via this exact same "type col row" format -
            // reuse it directly instead of re-reading the file here.
            if (!det.masked_pixels.empty()) {
                existing_masks_[name] = det.masked_pixels;
                WR_LOG(INFO, "Loaded " + std::to_string(existing_masks_[name].size()) + " existing masked pixels for " + name);
            }

            // Register plots, nested per-detector like the rest of the
            // framework, with names and binning matched to corry's MaskCreator.
            std::string dir = getName() + "/" + name;
            auto& pm = PlotManager::getInstance();

            pm.registerPlot2D(dir, "occupancy",
                det.n_pixels_x, -0.5, det.n_pixels_x - 0.5, det.n_pixels_y, -0.5, det.n_pixels_y - 0.5);
            pm.registerPlot2D(dir, "maskmap",
                det.n_pixels_x, -0.5, det.n_pixels_x - 0.5, det.n_pixels_y, -0.5, det.n_pixels_y - 0.5);

            if (method_ == "localdensity") {
                pm.registerPlot2D(dir, "density",
                    det.n_pixels_x, -0.5, det.n_pixels_x - 0.5, det.n_pixels_y, -0.5, det.n_pixels_y - 0.5);
                pm.registerPlot2D(dir, "local_significance",
                    det.n_pixels_x, -0.5, det.n_pixels_x - 0.5, det.n_pixels_y, -0.5, det.n_pixels_y - 0.5);
            }

            if (existing_masks_.count(name)) {
                Eigen::MatrixXd seed = Eigen::MatrixXd::Zero(det.n_pixels_x, det.n_pixels_y);
                for (auto const& p : existing_masks_[name]) seed(p.first, p.second) = 1.0;
                pm.getPlot2D(dir, "maskmap").setData(seed);
            }

            WR_LOG(DEBUG, "Masking enabled and plots registered for: " + name);
        }
    }

    void MaskCreator::run(DataBatch& batch) {
        total_events_ += global_config.get<size_t>("batch_size", 5000);

        std::vector<std::future<void>> futures;
        for (auto& entry : batch.pixels) {
            const std::string det_name = entry.first;
            auto& pixels = entry.second;

            // Skip if the detector type was filtered out during initialization
            if (occupancy_matrices_.find(det_name) == occupancy_matrices_.end()) continue;

            futures.push_back(thread_pool->submit([this, det_name, &pixels]() {
                auto& det = GeometryManager::getInstance().getDetector(det_name);

                // Accumulated into a thread-local matrix first (one task per
                // detector), merged into occupancy_matrices_ under that
                // detector's own lock only once at the end - avoids locking
                // per pixel.
                Eigen::MatrixXd local_occ = Eigen::MatrixXd::Zero(det.n_pixels_x, det.n_pixels_y);
                for (const auto& px : pixels) {
                    int c = px->column();
                    int r = px->row();
                    if (c >= 0 && c < det.n_pixels_x && r >= 0 && r < det.n_pixels_y) {
                        local_occ(c, r) += 1.0;
                    }
                }

                {
                    std::lock_guard<std::mutex> lock(*detector_mutexes_[det_name]);
                    occupancy_matrices_[det_name] += local_occ;
                }

                PlotManager::getInstance().getPlot2D(getName() + "/" + det_name, "occupancy").fillBatch(local_occ);
            }));
        }
        for (auto& f : futures) f.get();
    }

    void MaskCreator::finalize() {
        auto& pm = PlotManager::getInstance();

        for (auto& [name, matrix] : occupancy_matrices_) {
            auto& det = GeometryManager::getInstance().getDetector(name);
            std::string dir = getName() + "/" + name;

            Eigen::MatrixXd mask_display = Eigen::MatrixXd::Zero(det.n_pixels_x, det.n_pixels_y);
            std::vector<std::pair<int, int>> masked_indices;

            // Seed with pixels already masked via a pre-existing mask_file so
            // this run's output preserves them, matching corry.
            auto const& existing = existing_masks_[name];
            for (auto const& p : existing) {
                masked_indices.push_back(p);
                mask_display(p.first, p.second) = 1.0;
            }

            if (method_ == "frequency") {
                double mean = matrix.mean();
                for (int i = 0; i < matrix.rows(); ++i) {
                    for (int j = 0; j < matrix.cols(); ++j) {
                        if (existing.count({i, j})) continue;
                        if (matrix(i, j) > frequency_cut_ * mean) {
                            masked_indices.push_back({i, j});
                            mask_display(i, j) = 1.0;
                        }
                    }
                }
                WR_LOG(STATUS, name + ": Calculated Frequency mask (Cut: " + std::to_string(frequency_cut_) + " x Mean).");
            }
            else if (method_ == "localdensity") {
                // Kernel Density Estimation is performed in parallel sharded by rows.
                // Populates masked_indices/mask_display and the density/significance
                // diagnostic plots from the same computation used for the cut.
                calculateDensityMask(name, matrix, masked_indices, mask_display);
            }

            if (mask_dead_pixels_) {
                int dead_count = 0;
                for (int i = 0; i < matrix.rows(); ++i) {
                    for (int j = 0; j < matrix.cols(); ++j) {
                        if (matrix(i, j) == 0 && mask_display(i, j) == 0) {
                            mask_display(i, j) = 1.0;
                            masked_indices.push_back({i, j});
                            dead_count++;
                        }
                    }
                }
                WR_LOG(INFO, name + ": Identified " + std::to_string(dead_count) + " dead pixels.");
            }

            pm.getPlot2D(dir, "maskmap").setData(mask_display);

            std::filesystem::create_directories(getOutputPath());
            std::filesystem::path file_path = getOutputPath() / ("mask_" + name + ".txt");
            std::ofstream file(file_path);
            for (const auto& p : masked_indices) {
                file << "p " << p.first << " " << p.second << "\n";
            }

            WR_LOG(STATUS, name + ": Successfully masked " + std::to_string(masked_indices.size()) +
                           " total pixels. Output: " + file_path.string());
        }
    }

    void MaskCreator::calculateDensityMask(const std::string& det_name, const Eigen::MatrixXd& matrix,
                                            std::vector<std::pair<int, int>>& masked_indices, Eigen::MatrixXd& mask_display) {
        int rows = static_cast<int>(matrix.rows());
        int cols = static_cast<int>(matrix.cols());
        int bw = density_bandwidth_;

        Eigen::MatrixXd density = Eigen::MatrixXd::Zero(rows, cols);
        Eigen::MatrixXd significance = Eigen::MatrixXd::Zero(rows, cols);

        auto const& existing = existing_masks_[det_name];
        std::vector<std::pair<int, int>> hot_pixels;
        std::mutex hot_mtx;

        int n_threads = static_cast<int>(thread_pool->getThreadCount());
        int chunk = (rows + n_threads - 1) / n_threads;

        std::vector<std::future<void>> futures;
        for (int r_start = 0; r_start < rows; r_start += chunk) {
            int r_end = std::min(r_start + chunk, rows);

            // Each task only writes rows [r_start, r_end) of density/significance,
            // so no locking is needed for those - only the shared hot_pixels list.
            futures.push_back(thread_pool->submit([&, r_start, r_end, rows, cols, bw]() {
                std::vector<std::pair<int, int>> local_hot;
                for (int i = r_start; i < r_end; ++i) {
                    for (int j = 0; j < cols; ++j) {
                        double val = matrix(i, j);
                        if (val == 0) continue;

                        int i_min = std::max(0, i - bw), i_max = std::min(rows - 1, i + bw);
                        int j_min = std::max(0, j - bw), j_max = std::min(cols - 1, j + bw);

                        auto block = matrix.block(i_min, j_min, i_max - i_min + 1, j_max - j_min + 1);
                        double local_avg = (block.sum() - val) / (block.size() - 1);

                        double sig = (local_avg > 0) ? (val - local_avg) / std::sqrt(local_avg) : 0.0;
                        double rate = val / static_cast<double>(total_events_);

                        density(i, j) = local_avg;
                        significance(i, j) = sig;

                        if (!existing.count({i, j}) && (sig > sigma_max_ || rate > rate_max_)) {
                            local_hot.push_back({i, j});
                        }
                    }
                }
                std::lock_guard<std::mutex> lock(hot_mtx);
                hot_pixels.insert(hot_pixels.end(), local_hot.begin(), local_hot.end());
            }));
        }
        for (auto& f : futures) f.get();

        for (auto& p : hot_pixels) {
            masked_indices.push_back(p);
            mask_display(p.first, p.second) = 1.0;
        }

        auto& pm = PlotManager::getInstance();
        std::string dir = getName() + "/" + det_name;
        pm.getPlot2D(dir, "density").setData(density);
        pm.getPlot2D(dir, "local_significance").setData(significance);

        // Matches corry: for the local-density method the saved occupancy plot
        // is rescaled to a per-event rate rather than left as raw hit counts.
        if (total_events_ > 0) {
            pm.getPlot2D(dir, "occupancy").setData(matrix / static_cast<double>(total_events_));
        }

        WR_LOG(STATUS, det_name + ": Calculated Local Density mask (" + std::to_string(hot_pixels.size()) + " new pixels).");
    }

    REGISTER_MODULE(MaskCreator)
}