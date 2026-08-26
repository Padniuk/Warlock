/** @file AlignmentCAEN.cpp */
#include "AlignmentCAEN.hpp"
#include "core/ModuleFactory.hpp"
#include "core/detector/GeometryManager.hpp"
#include "core/utils/Logger.hpp"
#include "core/utils/PlotManager.hpp"
#include <algorithm>
#include <cmath>
#include <future>
#include <limits>
#include <set>

namespace framework {

    namespace {
        /// Zero-pads `n` to `width` digits (e.g. zeroPad(3, 2) -> "03") -
        /// keeps per-iteration plot directory names (iter01, iter02, ...,
        /// iter12) sorting in the right order instead of lexicographically
        /// (iter1, iter10, iter11, ...).
        std::string zeroPad(int n, int width) {
            std::string s = std::to_string(n);
            return s.size() < static_cast<size_t>(width) ? std::string(static_cast<size_t>(width) - s.size(), '0') + s
                                                           : s;
        }

        /// Maps `value` to a 0-based bin index over [lo, hi) split into
        /// `n_bins` equal-width bins, or -1 if out of range.
        int binIndex(double value, double lo, double hi, int n_bins) {
            if (n_bins <= 0 || hi <= lo) return -1;
            int idx = static_cast<int>((value - lo) / (hi - lo) * n_bins);
            if (idx < 0 || idx >= n_bins) return -1;
            return idx;
        }

        /// Weighted centroid of every bin of `passed` falling within one
        /// nominal pixel's local-frame footprint (`pixel_x`/`pixel_y` +/-
        /// half the pitch on each axis), weighted by each bin's own
        /// content. Returns false if the pixel's footprint doesn't
        /// overlap the map at all, or carries too little total weight
        /// (`min_weight`) to trust. Same algorithm as AnalysisEfficiency.cpp's
        /// identically-named helper.
        bool pixelCentroid(const Eigen::MatrixXd& passed, double x_lo, double x_hi, double y_lo, double y_hi,
                            double pixel_x, double pixel_y, double pitch_x, double pitch_y,
                            double min_weight, Eigen::Vector2d* centroid_out) {
            int n_bins_x = static_cast<int>(passed.rows());
            int n_bins_y = static_cast<int>(passed.cols());
            if (n_bins_x <= 0 || n_bins_y <= 0) return false;
            double bin_w = (x_hi - x_lo) / n_bins_x;
            double bin_h = (y_hi - y_lo) / n_bins_y;

            auto clampBin = [](int i, int n) { return std::max(0, std::min(n - 1, i)); };
            int ix_lo = clampBin(static_cast<int>((pixel_x - 0.5 * pitch_x - x_lo) / bin_w), n_bins_x);
            int ix_hi = clampBin(static_cast<int>((pixel_x + 0.5 * pitch_x - x_lo) / bin_w), n_bins_x);
            int iy_lo = clampBin(static_cast<int>((pixel_y - 0.5 * pitch_y - y_lo) / bin_h), n_bins_y);
            int iy_hi = clampBin(static_cast<int>((pixel_y + 0.5 * pitch_y - y_lo) / bin_h), n_bins_y);

            double sw = 0, swx = 0, swy = 0;
            for (int ix = ix_lo; ix <= ix_hi; ++ix) {
                double bx = x_lo + (ix + 0.5) * bin_w;
                for (int iy = iy_lo; iy <= iy_hi; ++iy) {
                    double w = passed(ix, iy);
                    if (w <= 0.0) continue;
                    double by = y_lo + (iy + 0.5) * bin_h;
                    sw += w; swx += w * bx; swy += w * by;
                }
            }
            if (sw < min_weight) return false;
            *centroid_out = {swx / sw, swy / sw};
            return true;
        }

        /// Closed-form 2D Procrustes rotation: the least-squares-optimal
        /// angle mapping `nominal` onto `measured` (both already de-meaned
        /// by the caller, so only rotation - not translation - is fit).
        /// theta = atan2(sum of 2D cross products, sum of dot products).
        double fitRotation(const std::vector<Eigen::Vector2d>& nominal, const std::vector<Eigen::Vector2d>& measured) {
            double cross_sum = 0.0, dot_sum = 0.0;
            for (size_t i = 0; i < nominal.size(); ++i) {
                cross_sum += nominal[i].x() * measured[i].y() - nominal[i].y() * measured[i].x();
                dot_sum += nominal[i].x() * measured[i].x() + nominal[i].y() * measured[i].y();
            }
            return std::atan2(cross_sum, dot_sum);
        }
    }

    AlignmentCAEN::AlignmentCAEN(Configuration cfg, Configuration g_cfg, ThreadPool* pool)
        : Module(std::move(cfg), std::move(g_cfg), pool) {
        dut_name_ = config.get<std::string>("name", "");
        type_filter_ = config.get<std::string>("type", "");

        // Cut-flow parameters, matching AnalysisEfficiency's own config
        // keys/defaults exactly, so a config block can be moved here
        // without changing which tracks count.
        chi2_cut_ = config.get<double>("chi2ndof_cut", 3.0);
        required_detectors_ = config.getArray<std::string>("require_associated_cluster_on");
        masked_pixel_distance_cut_ = config.get<int>("masked_pixel_distance_cut", 1);
        spatial_cut_sensoredge_ = config.get<double>("spatial_cut_sensoredge", 1.0);

        update_xy_translation_ = config.get<bool>("update_xy_translation", false);
        update_z_orientation_ = config.get<bool>("update_z_orientation", false);
        z_rotation_iterations_ = config.get<int>("z_rotation_iterations", 12);
        z_rotation_row_ = config.get<int>("z_rotation_row", 0);
        min_quadrant_entries_ = config.get<double>("min_quadrant_entries", 50.0);
        max_pixel_centroid_deviation_ = config.get<double>("max_pixel_centroid_deviation", 0.5);

        auto bin_sizes = config.getArray<double>("efficiency_2D_histograms_bin_size");
        bin_size_x_ = bin_sizes.size() >= 1 ? bin_sizes[0] : 0.01;
        bin_size_y_ = bin_sizes.size() >= 2 ? bin_sizes[1] : bin_size_x_;
    }

    void AlignmentCAEN::initialize() {
        auto& geo = GeometryManager::getInstance();

        // Same "name" (exactly one detector) / "type" (every matching
        // detector) resolution as DUTAssociation/AnalysisEfficiency.
        if (!dut_name_.empty()) {
            if (!geo.hasDetector(dut_name_)) {
                throw std::runtime_error("AlignmentCAEN: detector '" + dut_name_ + "' not found in geometry.");
            }
            target_detectors_.push_back(dut_name_);
        } else if (!type_filter_.empty()) {
            for (const auto& name : geo.getDetectorNames()) {
                if (geo.getDetector(name).type == type_filter_) target_detectors_.push_back(name);
            }
        } else {
            throw std::runtime_error("AlignmentCAEN: must configure either 'name' or 'type'.");
        }

        auto& pm = PlotManager::getInstance();
        for (const auto& det_name : target_detectors_) {
            auto& det = geo.getDetector(det_name);
            WR_LOG(STATUS, "Initializing CAEN alignment for DUT: " + det_name);

            DetAccum acc;
            acc.size_x = det.sizeX();
            acc.size_y = det.sizeY();
            acc.n_bins_x = std::max(1, static_cast<int>(std::round(1.4 * acc.size_x / bin_size_x_)));
            acc.n_bins_y = std::max(1, static_cast<int>(std::round(1.4 * acc.size_y / bin_size_y_)));
            acc.eff_num_local = Eigen::MatrixXd::Zero(acc.n_bins_x, acc.n_bins_y);
            accum_[det_name] = std::move(acc);

            // Diagnostic plots for the alignment principle itself, not the
            // fit result's own statistics (finalize()'s WR_LOG lines report
            // the actual corrections applied) - filled in finalize().
            std::string dir = getName() + "/" + det_name;
            if (update_xy_translation_) {
                double abs_sx = std::abs(det.sizeX()), abs_sy = std::abs(det.sizeY());
                // Matched-track global intercepts pooled into the XY mean -
                // this is the translation fit's entire input; the
                // histogram's own Mean (once converted to ROOT) is the
                // number written into det.position.
                pm.registerPlot1D(dir, "xy_translation_matched_x", 200,
                                   det.position.x() - abs_sx, det.position.x() + abs_sx);
                pm.registerPlot1D(dir, "xy_translation_matched_y", 200,
                                   det.position.y() - abs_sy, det.position.y() + abs_sy);
                pm.registerPlot2D(dir, "xy_translation_matched_2D", 200,
                                   det.position.x() - abs_sx, det.position.x() + abs_sx, 200,
                                   det.position.y() - abs_sy, det.position.y() + abs_sy);
            }
            if (update_z_orientation_) {
                // One point per trusted physical pixel, from the final
                // (converged) iteration: "nominal" is the pixel's position
                // on the DUT's own as-designed grid, "measured" is the
                // weighted centroid of matched tracks - overlaying the two
                // graphs shows the tilt fitRotation() solved for directly.
                // "convergence" plots the cumulative correction (degrees)
                // per iteration; a diverging fit (see the WARNING path
                // below) visibly fails to flatten out.
                pm.registerGraph(dir, "z_rotation_nominal_pixels");
                pm.registerGraph(dir, "z_rotation_measured_pixels");
                pm.registerGraph(dir, "z_rotation_convergence");
                // The actual data each pixel's "measured" point is a
                // weighted average of: matched-track local intercepts,
                // binned exactly like AnalysisEfficiency's own
                // localEfficiencyMap_trackPos. Lets the nominal/measured
                // pixel graphs above be checked directly against the
                // density map they were computed from, instead of taking
                // pixelCentroid()'s averaging on faith.
                double abs_sx = std::abs(det.sizeX()), abs_sy = std::abs(det.sizeY());
                const auto& acc_ref = accum_[det_name];
                pm.registerPlot2D(dir, "z_rotation_track_density_local", acc_ref.n_bins_x, -0.7 * abs_sx, 0.7 * abs_sx,
                                   acc_ref.n_bins_y, -0.7 * abs_sy, 0.7 * abs_sy);

                // Same three plots as above, once per iteration attempt
                // (not just the final one) - registered up front for the
                // full configured iteration cap, since finalize() doesn't
                // know in advance how many attempts a given DUT will
                // actually take before converging or being stopped as
                // diverging. Lets the fit's progress be watched step by
                // step instead of only seeing the end state.
                for (int iter = 1; iter <= z_rotation_iterations_; ++iter) {
                    std::string iter_dir = dir + "/z_rotation_steps/iter" + zeroPad(iter, 2);
                    pm.registerGraph(iter_dir, "nominal_pixels");
                    pm.registerGraph(iter_dir, "measured_pixels");
                    pm.registerPlot2D(iter_dir, "track_density_local", acc_ref.n_bins_x, -0.7 * abs_sx, 0.7 * abs_sx,
                                       acc_ref.n_bins_y, -0.7 * abs_sy, 0.7 * abs_sy);
                }
            }
        }
    }

    void AlignmentCAEN::run(DataBatch& batch) {
        if (batch.tracks.empty()) return;
        // Nothing to fit - don't bother accumulating records no one will read.
        if (!update_xy_translation_ && !update_z_orientation_) return;

        for (const auto& det_name : target_detectors_) {
            runTrackLoop(batch, det_name);
        }
    }

    void AlignmentCAEN::runTrackLoop(DataBatch& batch, const std::string& det_name) {
        auto& geo = GeometryManager::getInstance();
        auto& det = geo.getDetector(det_name);

        size_t n_tracks = batch.tracks.size();
        size_t n_threads = thread_pool->getThreadCount();
        size_t chunk = (n_tracks + n_threads - 1) / n_threads;
        std::vector<std::future<void>> futures;

        for (size_t i = 0; i < n_tracks; i += chunk) {
            size_t end = std::min(i + chunk, n_tracks);
            futures.push_back(thread_pool->submit([&, i, end]() {
                std::vector<TrackRecord> loc_records;

                for (size_t t_idx = i; t_idx < end; ++t_idx) {
                    auto const& track = batch.tracks[t_idx];

                    Eigen::Vector3d global_pred = track->interceptAt(det);
                    if (std::isnan(global_pred.x()) || std::isnan(global_pred.y())) continue;
                    Eigen::Vector3d local_pred = det.R_inv * (global_pred - det.position);

                    if (!det.hasIntercept(local_pred.x(), local_pred.y(), spatial_cut_sensoredge_)) continue;
                    if (track->getChi2ndof() > chi2_cut_) continue;
                    if (!track->isIsolated()) continue;
                    if (det.hitMasked(local_pred.x(), local_pred.y(), masked_pixel_distance_cut_)) continue;

                    bool requirements_met = true;
                    for (const auto& req : required_detectors_) {
                        bool has_req = false;
                        for (auto const& c : track->getAssociatedClusters()) {
                            if (c->detectorID() == req) { has_req = true; break; }
                        }
                        if (!has_req) { requirements_met = false; break; }
                    }
                    if (!requirements_met) continue;

                    bool has_associated = false;
                    for (auto const& c : track->getAssociatedClusters()) {
                        if (c->detectorID() == det_name) { has_associated = true; break; }
                    }

                    loc_records.push_back({global_pred, has_associated});
                }

                std::lock_guard<std::mutex> lock(acc_mtx_);
                auto& shared = accum_.at(det_name);
                shared.rotation_records.insert(shared.rotation_records.end(), loc_records.begin(), loc_records.end());
            }));
        }
        for (auto& f : futures) f.get();
    }

    void AlignmentCAEN::finalize() {
        if (!update_xy_translation_ && !update_z_orientation_) return;

        auto& geo = GeometryManager::getInstance();

        for (const auto& det_name : target_detectors_) {
            auto& acc = accum_.at(det_name);
            auto& det = geo.getDetector(det_name);

            double abs_sx = std::abs(acc.size_x), abs_sy = std::abs(acc.size_y);
            double x_lo = -0.7 * abs_sx, x_hi = 0.7 * abs_sx, y_lo = -0.7 * abs_sy, y_hi = 0.7 * abs_sy;

            // Re-derives local_pred for every stored matched record from
            // the detector's CURRENT (possibly just-corrected) R_inv/
            // position, and re-bins eff_num_local from scratch - lets this
            // be called again after orientation/position changes.
            auto rebinLocal = [&]() {
                acc.eff_num_local.setZero();
                for (auto const& rec : acc.rotation_records) {
                    if (!rec.matched) continue;
                    Eigen::Vector3d local_pred = det.R_inv * (rec.global_pred - det.position);
                    int lcol = binIndex(local_pred.x(), x_lo, x_hi, acc.n_bins_x);
                    int lrow = binIndex(local_pred.y(), y_lo, y_hi, acc.n_bins_y);
                    if (lcol < 0 || lrow < 0) continue;
                    acc.eff_num_local(lcol, lrow) += 1;
                }
            };

            // Per-pixel trust+plausibility, computed once against the
            // ORIGINAL geometry (before translation moves det.position)
            // and reused for both corrections below. Applied per physical
            // pixel rather than per row, so a single-row board (no "row"
            // concept when hasDieSplit() is false) gets the same
            // protection as a two-row board. See #max_pixel_centroid_deviation_'s
            // docs for what this catches.
            rebinLocal();
            std::vector<std::vector<bool>> pixel_ok(static_cast<size_t>(det.n_pixels_x),
                                                     std::vector<bool>(static_cast<size_t>(det.n_pixels_y), true));
            int n_bad_pixels = 0;
            for (int col = 0; col < det.n_pixels_x; ++col) {
                for (int row = 0; row < det.n_pixels_y; ++row) {
                    double px = det.getLocalX(col), py = det.getLocalY(row);
                    Eigen::Vector2d centroid;
                    bool trusted = pixelCentroid(acc.eff_num_local, x_lo, x_hi, y_lo, y_hi, px, py,
                                                  det.pitch_x, det.pitch_y, min_quadrant_entries_, &centroid);
                    bool plausible = trusted && (centroid - Eigen::Vector2d(px, py)).norm() <=
                                                     max_pixel_centroid_deviation_ * std::min(det.pitch_x, det.pitch_y);
                    // z_rotation_row_ is honored as an explicit manual
                    // override on top of the automatic per-pixel check -
                    // sometimes a whole die is known bad even when its
                    // statistics happen to look locally plausible.
                    bool row_allowed = true;
                    if (z_rotation_row_ > 0) row_allowed = (row == det.n_pixels_y - 1);
                    else if (z_rotation_row_ < 0) row_allowed = (row == 0);

                    bool ok = plausible && row_allowed;
                    pixel_ok[static_cast<size_t>(col)][static_cast<size_t>(row)] = ok;
                    if (!ok) {
                        ++n_bad_pixels;
                        WR_LOG(WARNING, det_name + ": pixel (col=" + std::to_string(col) + ", row=" + std::to_string(row) +
                                             ") excluded from alignment - " +
                                             (!row_allowed ? "excluded by z_rotation_row override"
                                              : trusted     ? "centroid implausibly far from nominal position"
                                                             : "too few matched tracks to trust") +
                                             ".");
                    }
                }
            }

            // For a die-split board, one bad pixel disqualifies its whole
            // die from the fit, not just itself: the two dies are
            // physically separate sensors with no guarantee they share a
            // single true tilt, so mixing a bad pixel's healthy sibling(s)
            // into the same rotation fit as the other die would violate
            // that guarantee. The per-pixel check above still runs first,
            // since a single-row board has no other row to fall back on.
            if (det.hasDieSplit()) {
                std::set<std::string> bad_dies;
                for (int col = 0; col < det.n_pixels_x; ++col) {
                    for (int row = 0; row < det.n_pixels_y; ++row) {
                        if (!pixel_ok[static_cast<size_t>(col)][static_cast<size_t>(row)]) {
                            std::string die = det.dieOf(col, row);
                            if (!die.empty()) bad_dies.insert(die);
                        }
                    }
                }
                for (int col = 0; col < det.n_pixels_x; ++col) {
                    for (int row = 0; row < det.n_pixels_y; ++row) {
                        if (!pixel_ok[static_cast<size_t>(col)][static_cast<size_t>(row)]) continue;
                        std::string die = det.dieOf(col, row);
                        if (bad_dies.count(die)) {
                            pixel_ok[static_cast<size_t>(col)][static_cast<size_t>(row)] = false;
                            ++n_bad_pixels;
                            WR_LOG(WARNING, det_name + ": pixel (col=" + std::to_string(col) + ", row=" + std::to_string(row) +
                                                 ") excluded from alignment - shares die '" + die +
                                                 "' with an already-excluded pixel, and this board's dies aren't "
                                                 "guaranteed to share a single true tilt.");
                        }
                    }
                }
            }

            // If exactly one die's pixels survived #pixel_ok (the other
            // excluded entirely, whether by die-propagation above or by an
            // explicit z_rotation_row_ override), the translation mean
            // below is computed only from that die's own matched tracks,
            // so it lands on the surviving die's own center rather than
            // the whole board's reference origin. Records which row is the
            // surviving die so the translation block can shift the raw
            // mean back by that die's own nominal offset (getLocalY())
            // before writing det.position. -1 means no single-die
            // compensation is needed (a non-split board, or a split board
            // where both dies still have at least one surviving pixel).
            int healthy_row_for_translation = -1;
            if (det.hasDieSplit() && det.n_pixels_y >= 2) {
                bool bottom_ok = false, top_ok = false;
                for (int col = 0; col < det.n_pixels_x; ++col) {
                    if (pixel_ok[static_cast<size_t>(col)][0]) bottom_ok = true;
                    if (pixel_ok[static_cast<size_t>(col)][static_cast<size_t>(det.n_pixels_y - 1)]) top_ok = true;
                }
                if (bottom_ok && !top_ok) healthy_row_for_translation = 0;
                else if (!bottom_ok && top_ok) healthy_row_for_translation = det.n_pixels_y - 1;
            }

            // Nearest physical pixel to a local-frame position, clamped to
            // the board's actual pixel range - used to look up #pixel_ok
            // for a matched track's own intercept. The clamp is a safety
            // fallback; a track passing DUTAssociation's cut to get here
            // is already close to a real pixel.
            auto pixelOfLocal = [&](double lx, double ly) -> std::pair<int, int> {
                int col = static_cast<int>(std::lround(det.getColumn(lx)));
                int row = static_cast<int>(std::lround(det.getRow(ly)));
                col = std::max(0, std::min(det.n_pixels_x - 1, col));
                row = std::max(0, std::min(det.n_pixels_y - 1, row));
                return {col, row};
            };

            // Target position is the mean of matched global intercepts in
            // LOCAL coordinates - not the track-vs-cluster residual: that
            // residual, R_inv*(global_pred - global_cluster), cancels the
            // position terms exactly, so it's independent of det.position
            // and can't drive a translation correction. Single pass, not
            // iterative: this mean doesn't depend on det.position, so
            // recomputing it after moving position gives the same answer.
            if (update_xy_translation_) {
                auto& pm = PlotManager::getInstance();
                std::string dir = getName() + "/" + det_name;
                double sum_x = 0.0, sum_y = 0.0;
                size_t count = 0, excluded = 0;
                for (auto const& rec : acc.rotation_records) {
                    if (!rec.matched) continue;
                    Eigen::Vector3d local_pred = det.R_inv * (rec.global_pred - det.position);
                    auto [col, row] = pixelOfLocal(local_pred.x(), local_pred.y());
                    if (!pixel_ok[static_cast<size_t>(col)][static_cast<size_t>(row)]) { ++excluded; continue; }
                    sum_x += rec.global_pred.x();
                    sum_y += rec.global_pred.y();
                    ++count;
                    // The translation fit's actual input, point by point -
                    // these histograms' own Mean is literally the number
                    // written into det.position below.
                    pm.fill1D(dir, "xy_translation_matched_x", rec.global_pred.x());
                    pm.fill1D(dir, "xy_translation_matched_y", rec.global_pred.y());
                    pm.fill2D(dir, "xy_translation_matched_2D", rec.global_pred.x(), rec.global_pred.y());
                }

                if (count < 100) {
                    WR_LOG(WARNING, det_name + ": only " + std::to_string(count) +
                                         " matched track(s) after excluding " + std::to_string(excluded) +
                                         " landing on untrusted pixels - too few to fit a translation, left uncorrected.");
                } else {
                    double mean_x = sum_x / static_cast<double>(count);
                    double mean_y = sum_y / static_cast<double>(count);
                    if (healthy_row_for_translation >= 0) {
                        // mean_x/mean_y above is the surviving die's own
                        // centroid, not the whole board's - shift back by
                        // that die's nominal offset from the board's local
                        // origin (getLocalY()) so the healthy die lands
                        // where its own data says, leaving the excluded
                        // die at its fixed, configured relative spacing.
                        Eigen::Vector3d nominal_local(0.0, det.getLocalY(healthy_row_for_translation), 0.0);
                        Eigen::Vector3d origin = Eigen::Vector3d(mean_x, mean_y, det.position.z()) - det.R * nominal_local;
                        mean_x = origin.x();
                        mean_y = origin.y();
                    }
                    WR_LOG(STATUS, det_name + " translation aligned -> center: (" + std::to_string(mean_x) +
                                        ", " + std::to_string(mean_y) + ")mm, from " + std::to_string(count) +
                                        " matched track(s), " + std::to_string(excluded) +
                                        " excluded from untrusted pixels" +
                                        (healthy_row_for_translation >= 0 ? " (healthy die only, offset-compensated)" : "") +
                                        " [was: (" + std::to_string(det.position(0)) +
                                        ", " + std::to_string(det.position(1)) + ")mm]");
                    det.position(0) = mean_x;
                    det.position(1) = mean_y;
                }
                // No updateMatrices() call needed here - R/R_inv/V_inv are
                // derived from orientation and spatial_res_x/y only, never
                // from position, so a pure translation change doesn't
                // invalidate them.
            }

            // One (nominal, measured) point pair per physical pixel that
            // #pixel_ok trusts - nominal is the pixel's position on the
            // DUT's own (unrotated) grid, measured is the weighted
            // centroid of matched tracks in its footprint, both in this
            // detector's local frame. Iterated: a rotation changes which
            // bin (and pixel) each stored intercept falls into, so measured
            // centroids need re-binning against the corrected geometry
            // each pass - unlike translation, this fit depends on the
            // current det.orientation.
            if (update_z_orientation_) {
                auto& pm = PlotManager::getInstance();
                std::string dir = getName() + "/" + det_name;
                // Raw (pre-demean) point pairs, overwritten every
                // iteration - after the loop this holds the final,
                // converged (or last-attempted) state, used for the
                // nominal/measured graphs below regardless of which break
                // path ended the loop.
                std::vector<Eigen::Vector2d> last_nominal_raw, last_measured_raw;
                double last_dtheta = 0.0; ///< The specific correction that produced last_nominal_raw/last_measured_raw, for the residual plots below.

                double total_dtheta = 0.0;
                int iterations_run = 0;
                int iterations_attempted = 0; ///< Snapshotted attempts (see z_rotation_steps/ below), including one that ends up rejected by the divergence guard - used after the loop to discard the unused tail of pre-registered iteration slots.
                double prev_abs_dtheta = std::numeric_limits<double>::infinity();
                for (int iter = 0; iter < z_rotation_iterations_; ++iter) {
                    rebinLocal();

                    std::vector<Eigen::Vector2d> nominal, measured;
                    for (int col = 0; col < det.n_pixels_x; ++col) {
                        for (int row = 0; row < det.n_pixels_y; ++row) {
                            if (!pixel_ok[static_cast<size_t>(col)][static_cast<size_t>(row)]) continue;
                            double px = det.getLocalX(col), py = det.getLocalY(row);
                            Eigen::Vector2d centroid;
                            if (!pixelCentroid(acc.eff_num_local, x_lo, x_hi, y_lo, y_hi, px, py,
                                                det.pitch_x, det.pitch_y, min_quadrant_entries_, &centroid)) {
                                continue;
                            }
                            nominal.push_back({px, py});
                            measured.push_back(centroid);
                        }
                    }

                    if (nominal.size() < 2) {
                        WR_LOG(WARNING, det_name + ": only " + std::to_string(nominal.size()) +
                                             " pixel(s) with enough passed-histogram statistics (need >=2) - " +
                                             "stopped after " + std::to_string(iterations_run) + " iteration(s).");
                        break;
                    }

                    // Saved before de-meaning mutates nominal/measured in
                    // place below - only committed to last_nominal_raw/
                    // last_measured_raw once this iteration's correction is
                    // actually applied, so they stay in sync with whatever
                    // geometry is kept, including when the divergence guard
                    // rejects the last attempted iteration.
                    std::vector<Eigen::Vector2d> raw_nominal = nominal, raw_measured = measured;

                    // Snapshot this attempt (whether or not it ends up
                    // accepted below) into its own per-iteration plots -
                    // acc.eff_num_local still holds the density map
                    // rebinLocal() just rebuilt at the top of this pass.
                    {
                        std::string iter_dir = dir + "/z_rotation_steps/iter" + zeroPad(iter + 1, 2);
                        for (auto const& p : raw_nominal) pm.addPoint(iter_dir, "nominal_pixels", p.x(), p.y());
                        for (auto const& p : raw_measured) pm.addPoint(iter_dir, "measured_pixels", p.x(), p.y());
                        pm.getPlot2D(iter_dir, "track_density_local").setData(acc.eff_num_local);
                        iterations_attempted = iter + 1;
                    }

                    // De-mean both point sets - removes any translation
                    // offset between the nominal grid and the measured
                    // centroids, leaving pure rotation for fitRotation().
                    Eigen::Vector2d nominal_mean = Eigen::Vector2d::Zero(), measured_mean = Eigen::Vector2d::Zero();
                    for (size_t i = 0; i < nominal.size(); ++i) { nominal_mean += nominal[i]; measured_mean += measured[i]; }
                    nominal_mean /= static_cast<double>(nominal.size());
                    measured_mean /= static_cast<double>(measured.size());
                    for (size_t i = 0; i < nominal.size(); ++i) { nominal[i] -= nominal_mean; measured[i] -= measured_mean; }

                    double dtheta = fitRotation(nominal, measured);

                    // Once a correction has settled to the noise floor
                    // (bin-quantization/floating-point scale, not a real
                    // remaining tilt), later iterations can tick up by a
                    // sub-microdegree amount purely from that noise - not
                    // genuine divergence. Apply this negligible correction
                    // and stop cleanly rather than letting the (much
                    // coarser) divergence check below misfire on it.
                    constexpr double kConvergedEpsilonRad = 1e-6;
                    if (std::abs(dtheta) < kConvergedEpsilonRad) {
                        det.orientation(2) += dtheta;
                        det.updateMatrices();
                        total_dtheta += dtheta;
                        ++iterations_run;
                        last_nominal_raw = raw_nominal;
                        last_measured_raw = raw_measured;
                        last_dtheta = dtheta;
                        pm.addPoint(dir, "z_rotation_convergence", iterations_run, total_dtheta * 180.0 / M_PI);
                        WR_LOG(DEBUG, det_name + ": Z-rotation converged after " + std::to_string(iterations_run) + " iteration(s).");
                        break;
                    }

                    // A converging fit's correction shrinks every
                    // iteration. Stopping as soon as a correction stops
                    // shrinking, rather than blindly running the full
                    // iteration count, keeps a diverging DUT's bad data
                    // from applying an ever-growing bogus rotation.
                    if (iter > 0 && std::abs(dtheta) >= prev_abs_dtheta) {
                        WR_LOG(WARNING, det_name + ": Z-rotation fit is diverging (iteration " + std::to_string(iterations_run + 1) +
                                             " correction " + std::to_string(dtheta * 180.0 / M_PI) + "deg is not smaller than the previous " +
                                             std::to_string(prev_abs_dtheta * 180.0 / M_PI) +
                                             "deg) - stopping early and keeping the result from " + std::to_string(iterations_run) +
                                             " prior iteration(s). This DUT's matched-track statistics may still be too poor/asymmetric "
                                             "for this fit to trust, even after excluding untrusted pixels.");
                        break;
                    }
                    prev_abs_dtheta = std::abs(dtheta);

                    det.orientation(2) += dtheta;
                    det.updateMatrices();
                    total_dtheta += dtheta;
                    ++iterations_run;
                    last_nominal_raw = raw_nominal;
                    last_measured_raw = raw_measured;
                    last_dtheta = dtheta;
                    pm.addPoint(dir, "z_rotation_convergence", iterations_run, total_dtheta * 180.0 / M_PI);
                    WR_LOG(DEBUG, det_name + ": Z-rotation iteration " + std::to_string(iterations_run) + "/" +
                                       std::to_string(z_rotation_iterations_) + " = " + std::to_string(dtheta * 180.0 / M_PI) +
                                       "deg, from " + std::to_string(nominal.size()) + " pixel centroid(s).");
                }

                // Discard the unused tail of pre-registered iteration
                // slots (z_rotation_iterations_ is a safe upper bound
                // fixed at initialize()-time, before how many attempts
                // this DUT would actually need was known) - otherwise
                // they'd still show up in the output as empty groups.
                for (int iter = iterations_attempted + 1; iter <= z_rotation_iterations_; ++iter) {
                    std::string iter_dir = dir + "/z_rotation_steps/iter" + zeroPad(iter, 2);
                    pm.unregisterGraph(iter_dir, "nominal_pixels");
                    pm.unregisterGraph(iter_dir, "measured_pixels");
                    pm.unregisterPlot2D(iter_dir, "track_density_local");
                }

                // One point per trusted pixel, from the final iteration -
                // overlaying these two graphs shows the tilt the fit
                // solved for directly, rather than just the resulting
                // angle.
                for (auto const& p : last_nominal_raw) pm.addPoint(dir, "z_rotation_nominal_pixels", p.x(), p.y());
                for (auto const& p : last_measured_raw) pm.addPoint(dir, "z_rotation_measured_pixels", p.x(), p.y());

                // acc.eff_num_local still holds the last iteration's
                // rebinLocal() state - the exact density map
                // pixelCentroid() averaged over to produce every
                // "measured" point above.
                pm.getPlot2D(dir, "z_rotation_track_density_local").setData(acc.eff_num_local);

                // Per-pixel agreement check for the fit that produced
                // last_dtheta: each trusted pixel implies its own rotation
                // angle from itself alone; how far that is from the actual
                // fitted angle, averaged over pixels, says how much the
                // pixels agree with each other - reported in the log line
                // below rather than as a plot (as few as 4 physical pixels
                // on a 2x2 board is too few points for any histogram or
                // scatter graph to usefully show).
                double sum_angle_resid = 0.0, sum_angle_resid2 = 0.0;
                size_t n_resid = 0;
                if (last_nominal_raw.size() >= 2) {
                    Eigen::Vector2d nom_mean = Eigen::Vector2d::Zero(), meas_mean = Eigen::Vector2d::Zero();
                    for (size_t i = 0; i < last_nominal_raw.size(); ++i) {
                        nom_mean += last_nominal_raw[i];
                        meas_mean += last_measured_raw[i];
                    }
                    nom_mean /= static_cast<double>(last_nominal_raw.size());
                    meas_mean /= static_cast<double>(last_measured_raw.size());

                    for (size_t i = 0; i < last_nominal_raw.size(); ++i) {
                        Eigen::Vector2d nom = last_nominal_raw[i] - nom_mean;
                        Eigen::Vector2d meas = last_measured_raw[i] - meas_mean;
                        if (nom.norm() < 1e-9) continue;

                        double cross_i = nom.x() * meas.y() - nom.y() * meas.x();
                        double dot_i = nom.x() * meas.x() + nom.y() * meas.y();
                        double angle_i = std::atan2(cross_i, dot_i);
                        double angle_resid_deg = (angle_i - last_dtheta) * 180.0 / M_PI;
                        sum_angle_resid += angle_resid_deg;
                        sum_angle_resid2 += angle_resid_deg * angle_resid_deg;
                        ++n_resid;
                    }
                }

                std::string resid_summary;
                if (n_resid >= 2) {
                    double mean_resid = sum_angle_resid / static_cast<double>(n_resid);
                    double var_resid = sum_angle_resid2 / static_cast<double>(n_resid) - mean_resid * mean_resid;
                    double std_resid = std::sqrt(std::max(0.0, var_resid));
                    resid_summary = ", per-pixel angle agreement = " + std::to_string(std_resid) +
                                     "deg RMS (n=" + std::to_string(n_resid) + ")";
                }

                WR_LOG(STATUS, det_name + ": total Z-rotation correction = " + std::to_string(total_dtheta * 180.0 / M_PI) +
                                   "deg over " + std::to_string(iterations_run) + " iteration(s), " +
                                   std::to_string(n_bad_pixels) + " pixel(s) excluded from the fit" + resid_summary + ".");
            }
        }

        std::string out_geo = global_config.get<std::string>("detectors_file_updated", "z_rotated.geo");
        geo.saveGeometry(out_geo);
        WR_LOG(STATUS, "AlignmentCAEN: corrected geometry saved to " + out_geo);
    }

    REGISTER_MODULE(AlignmentCAEN)
}
