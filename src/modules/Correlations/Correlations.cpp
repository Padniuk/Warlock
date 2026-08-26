#include "Correlations.hpp"
#include "core/ModuleFactory.hpp"
#include "core/utils/Logger.hpp"
#include "core/detector/GeometryManager.hpp"
#include "core/utils/PlotManager.hpp"
#include <unordered_map>
#include <vector>
#include <future>
#include <cmath>

namespace framework {

    namespace {
        // True (unbinned) first/second moments and cross-moment of every
        // fill2DLocal() call, regardless of in-range/out-of-range (matching
        // ROOT TH2::Fill()'s own unconditional statistics update) - passed
        // to Histogram2D's exact-moment fillBatch() overload so
        // getEntries()-adjacent ROOT statistics (Mean/StdDev, via
        // TH2::PutStats() in H5ToRootConverter) match corryvreckan's native
        // TH2::Fill() precisely instead of only approximating it from bin
        // centers - see PlotManager.hpp's Histogram2D docs. Plain summed
        // fields, no locking needed: one instance lives per-thread inside
        // DetAccum, only summed together at merge time (already
        // single-threaded there).
        struct Moments2D {
            double sw = 0, sx = 0, sy = 0, sx2 = 0, sy2 = 0, sxy = 0;
        };
        Moments2D& operator+=(Moments2D& a, const Moments2D& b) {
            a.sw += b.sw; a.sx += b.sx; a.sy += b.sy; a.sx2 += b.sx2; a.sy2 += b.sy2; a.sxy += b.sxy;
            return a;
        }

        // Mirrors Histogram2D::fill()'s exact bin-edge/width formula, so a
        // matrix built up here and merged via fillBatch() lands in the same
        // bins the single-threaded fill() calls would have used.
        void fill2DLocal(Eigen::MatrixXd& mat, Moments2D& mom, double x, double y,
                                 double minx, double maxx, double miny, double maxy) {
            // Matches ROOT TH2::Fill()'s default (StatOverflows() off): an
            // out-of-range pair still lands in the matrix's (conceptual)
            // overflow region but is excluded from the moments used for
            // GetMean()/GetStdDev() - see Histogram1D::fill()'s equivalent
            // comment in PlotManager.cpp.
            if (x < minx || x >= maxx || y < miny || y >= maxy) return;
            mom.sw += 1.0; mom.sx += x; mom.sy += y; mom.sx2 += x * x; mom.sy2 += y * y; mom.sxy += x * y;
            double wx = (maxx - minx) / static_cast<double>(mat.rows());
            double wy = (maxy - miny) / static_cast<double>(mat.cols());
            size_t ix = static_cast<size_t>((x - minx) / wx);
            size_t iy = static_cast<size_t>((y - miny) / wy);
            mat(ix, iy) += 1.0;
        }

        // Per-thread accumulator for one non-reference detector. Filling
        // PlotManager's shared histograms directly from every thread would
        // serialize on that histogram's own mutex (correlation plots are
        // filled once per hit/cluster PAIR, so that's a lot of lock
        // acquisitions under contention); accumulating locally and merging
        // once per batch via fillBatch() takes that lock only a handful of
        // times total instead.
        struct DetAccum {
            std::vector<double> corrX, corrY, corrXY, corrYX;
            std::vector<double> corrTime, corrTimeInt, corrTimePx;
            Eigen::MatrixXd corrX2Dlocal, corrY2Dlocal;
            Eigen::MatrixXd corrColCol_px, corrRowRow_px, corrColRow_px, corrRowCol_px;
            Eigen::MatrixXd corrX2D, corrY2D, corrXY2D, corrYX2D;
            Moments2D momX2Dlocal, momY2Dlocal;
            Moments2D momColCol_px, momRowRow_px, momColRow_px, momRowCol_px;
            Moments2D momX2D, momY2D, momXY2D, momYX2D;
        };
    }

    Correlations::Correlations(Configuration cfg, Configuration g_cfg, ThreadPool* pool)
        : Module(std::move(cfg), std::move(g_cfg), pool) {
        time_cut_abs_ = config.get<double>("time_cut_abs", 690000.0);
        do_time_cut_ = config.get<bool>("do_time_cut", false);
        time_binning_ = config.get<double>("time_binning", 1.0);
    }

    void Correlations::initialize() {
        WR_LOG(STATUS, "Initializing Correlations...");
        auto& geo = GeometryManager::getInstance();
        std::string ref = geo.getReferenceName();

        int time_bins = static_cast<int>(2.0 * time_cut_abs_ / time_binning_);
        double time_min = -time_cut_abs_ - time_binning_ / 2.0;
        double time_max = time_cut_abs_ - time_binning_ / 2.0;

        for (const auto& name : geo.getDetectorNames()) {
            auto& det = geo.getDetector(name);
            auto& ref_det = geo.getDetector(ref);
            // corry computes and stores a real self-correlation for the
            // reference detector too (non-trivial whenever an event has
            // more than one cluster on the reference plane) - it never
            // special-cases `name == ref`, so this doesn't either.
            if (!det.isPixelDetector()) continue;

            std::string dir = getName() + "/" + name;
            auto& pm = PlotManager::getInstance();

            // Binning matched exactly to corryvreckan's Correlations module.
            pm.registerPlot2D(dir, "hitmap", det.n_pixels_x, -0.5, det.n_pixels_x - 0.5, det.n_pixels_y, -0.5, det.n_pixels_y - 0.5);
            pm.registerPlot2D(dir, "hitmap_clusters", det.n_pixels_x, -0.5, det.n_pixels_x - 0.5, det.n_pixels_y, -0.5, det.n_pixels_y - 0.5);
            pm.registerPlot1D(dir, "eventTimes", 3000000, -1e-5, 300 - 1e-5);

            pm.registerPlot1D(dir, "correlationX", 1000, -10.01, 9.99);
            pm.registerPlot1D(dir, "correlationY", 1000, -10.01, 9.99);
            pm.registerPlot1D(dir, "correlationXY", 1000, -10.01, 9.99);
            pm.registerPlot1D(dir, "correlationYX", 1000, -10.01, 9.99);
            pm.registerPlot1D(dir, "correlationTime", time_bins, time_min, time_max);
            pm.registerPlot1D(dir, "correlationTime_px", time_bins, time_min, time_max);
            pm.registerPlot1D(dir, "correlationTimeInt", 8000, -40005, 39995);

            pm.registerPlot2D(dir, "correlationX_2Dlocal", det.n_pixels_x, -0.5, det.n_pixels_x - 0.5, ref_det.n_pixels_x, -0.5, ref_det.n_pixels_x - 0.5);
            pm.registerPlot2D(dir, "correlationY_2Dlocal", det.n_pixels_y, -0.5, det.n_pixels_y - 0.5, ref_det.n_pixels_y, -0.5, ref_det.n_pixels_y - 0.5);

            pm.registerPlot2D(dir, "correlationColCol_px", det.n_pixels_x, -0.5, det.n_pixels_x - 0.5, ref_det.n_pixels_x, -0.5, ref_det.n_pixels_x - 0.5);
            pm.registerPlot2D(dir, "correlationRowRow_px", det.n_pixels_y, -0.5, det.n_pixels_y - 0.5, ref_det.n_pixels_y, -0.5, ref_det.n_pixels_y - 0.5);
            pm.registerPlot2D(dir, "correlationColRow_px", det.n_pixels_x, -0.5, det.n_pixels_x - 0.5, ref_det.n_pixels_y, -0.5, ref_det.n_pixels_y - 0.5);
            pm.registerPlot2D(dir, "correlationRowCol_px", det.n_pixels_y, -0.5, det.n_pixels_y - 0.5, ref_det.n_pixels_x, -0.5, ref_det.n_pixels_x - 0.5);

            pm.registerPlot2D(dir, "correlationX_2D", 100, -10.1, 9.9, 100, -10.1, 9.9);
            pm.registerPlot2D(dir, "correlationY_2D", 100, -10.1, 9.9, 100, -10.1, 9.9);
            pm.registerPlot2D(dir, "correlationXY_2D", 100, -10.1, 9.9, 100, -10.1, 9.9);
            pm.registerPlot2D(dir, "correlationYX_2D", 100, -10.1, 9.9, 100, -10.1, 9.9);
        }
    }

    void Correlations::run(DataBatch& batch) {
        auto& pm = PlotManager::getInstance();
        auto& geo = GeometryManager::getInstance();
        std::string ref_name = geo.getReferenceName();

        // Raw hitmap/eventTimes/hitmap_clusters are unconditional on a
        // reference being present, matching corry filling them at the very
        // start of every event, independent of its own reference-pairing
        // loop below. This part is O(n_hits), not O(n_hits x n_ref_hits)
        // like the correlation pairing, so it stays as simple direct fills.
        for (auto const& [det_name, pixels] : batch.pixels) {
            std::string dir = getName() + "/" + det_name;
            if (!pm.hasPlot1D(dir, "correlationX")) continue;
            for (auto const& p : pixels) {
                pm.fill2D(dir, "hitmap", p->column(), p->row());
                pm.fill1D(dir, "eventTimes", p->timestamp() * 1e-9);
            }
        }
        for (auto const& [det_name, clusters] : batch.clusters) {
            std::string dir = getName() + "/" + det_name;
            if (!pm.hasPlot1D(dir, "correlationX")) continue;
            for (auto const& c : clusters) {
                pm.fill2D(dir, "hitmap_clusters", c->column(), c->row());
            }
        }

        size_t n_threads = thread_pool ? thread_pool->getThreadCount() : 1;
        if (n_threads < 1) n_threads = 1;

        // --- Pixel-level ("_px") correlations ---
        if (batch.pixels.count(ref_name) > 0) {
            std::unordered_map<uint64_t, std::unordered_map<std::string, std::vector<std::shared_ptr<Pixel>>>> pixel_events;
            for (auto const& [det, pixels] : batch.pixels) {
                for (auto const& p : pixels) pixel_events[p->eventID()][det].push_back(p);
            }

            std::vector<decltype(pixel_events)::value_type*> event_ptrs;
            event_ptrs.reserve(pixel_events.size());
            for (auto& kv : pixel_events) event_ptrs.push_back(&kv);

            std::vector<std::unordered_map<std::string, DetAccum>> thread_accums(n_threads);
            size_t chunk = (event_ptrs.size() + n_threads - 1) / std::max<size_t>(n_threads, 1);
            std::vector<std::future<void>> futures;

            for (size_t t = 0; t < n_threads; ++t) {
                size_t start = t * chunk, end = std::min(start + chunk, event_ptrs.size());
                if (start >= end) continue;
                futures.push_back(thread_pool->submit([&, start, end, t]() {
                    auto& local = thread_accums[t];
                    for (size_t k = start; k < end; ++k) {
                        auto& dets = event_ptrs[k]->second;
                        if (dets.count(ref_name) == 0) continue;
                        auto const& ref_pixels = dets.at(ref_name);

                        for (auto const& [name, pixels] : dets) {
                            std::string dir = getName() + "/" + name;
                            if (!pm.hasPlot1D(dir, "correlationX")) continue;

                            auto& det = geo.getDetector(name);
                            auto& ref_det = geo.getDetector(ref_name);
                            auto& acc = local[name];
                            if (acc.corrColCol_px.size() == 0) {
                                acc.corrColCol_px = Eigen::MatrixXd::Zero(det.n_pixels_x, ref_det.n_pixels_x);
                                acc.corrRowRow_px = Eigen::MatrixXd::Zero(det.n_pixels_y, ref_det.n_pixels_y);
                                acc.corrColRow_px = Eigen::MatrixXd::Zero(det.n_pixels_x, ref_det.n_pixels_y);
                                acc.corrRowCol_px = Eigen::MatrixXd::Zero(det.n_pixels_y, ref_det.n_pixels_x);
                            }

                            for (auto const& p : pixels) {
                                for (auto const& rp : ref_pixels) {
                                    fill2DLocal(acc.corrColCol_px, acc.momColCol_px, p->column(), rp->column(), -0.5, det.n_pixels_x - 0.5, -0.5, ref_det.n_pixels_x - 0.5);
                                    fill2DLocal(acc.corrRowRow_px, acc.momRowRow_px, p->row(), rp->row(), -0.5, det.n_pixels_y - 0.5, -0.5, ref_det.n_pixels_y - 0.5);
                                    fill2DLocal(acc.corrColRow_px, acc.momColRow_px, p->column(), rp->row(), -0.5, det.n_pixels_x - 0.5, -0.5, ref_det.n_pixels_y - 0.5);
                                    fill2DLocal(acc.corrRowCol_px, acc.momRowCol_px, p->row(), rp->column(), -0.5, det.n_pixels_y - 0.5, -0.5, ref_det.n_pixels_x - 0.5);
                                    acc.corrTimePx.push_back(rp->timestamp() - p->timestamp());
                                }
                            }
                        }
                    }
                }));
            }
            for (auto& f : futures) f.get();

            // Merge: one lock acquisition per plot per detector, not per pair.
            std::unordered_map<std::string, DetAccum> merged;
            for (auto& local : thread_accums) {
                for (auto& [name, acc] : local) {
                    auto& m = merged[name];
                    if (m.corrColCol_px.size() == 0) {
                        m.corrColCol_px = acc.corrColCol_px; m.corrRowRow_px = acc.corrRowRow_px;
                        m.corrColRow_px = acc.corrColRow_px; m.corrRowCol_px = acc.corrRowCol_px;
                    } else {
                        m.corrColCol_px += acc.corrColCol_px; m.corrRowRow_px += acc.corrRowRow_px;
                        m.corrColRow_px += acc.corrColRow_px; m.corrRowCol_px += acc.corrRowCol_px;
                    }
                    m.momColCol_px += acc.momColCol_px; m.momRowRow_px += acc.momRowRow_px;
                    m.momColRow_px += acc.momColRow_px; m.momRowCol_px += acc.momRowCol_px;
                    m.corrTimePx.insert(m.corrTimePx.end(), acc.corrTimePx.begin(), acc.corrTimePx.end());
                }
            }
            for (auto& [name, acc] : merged) {
                std::string dir = getName() + "/" + name;
                pm.getPlot2D(dir, "correlationColCol_px").fillBatch(acc.corrColCol_px, acc.momColCol_px.sw, acc.momColCol_px.sx,
                                    acc.momColCol_px.sy, acc.momColCol_px.sx2, acc.momColCol_px.sy2, acc.momColCol_px.sxy);
                pm.getPlot2D(dir, "correlationRowRow_px").fillBatch(acc.corrRowRow_px, acc.momRowRow_px.sw, acc.momRowRow_px.sx,
                                    acc.momRowRow_px.sy, acc.momRowRow_px.sx2, acc.momRowRow_px.sy2, acc.momRowRow_px.sxy);
                pm.getPlot2D(dir, "correlationColRow_px").fillBatch(acc.corrColRow_px, acc.momColRow_px.sw, acc.momColRow_px.sx,
                                    acc.momColRow_px.sy, acc.momColRow_px.sx2, acc.momColRow_px.sy2, acc.momColRow_px.sxy);
                pm.getPlot2D(dir, "correlationRowCol_px").fillBatch(acc.corrRowCol_px, acc.momRowCol_px.sw, acc.momRowCol_px.sx,
                                    acc.momRowCol_px.sy, acc.momRowCol_px.sx2, acc.momRowCol_px.sy2, acc.momRowCol_px.sxy);
                pm.getPlot1D(dir, "correlationTime_px").fillBatch(acc.corrTimePx);
            }
        }

        if (batch.clusters.count(ref_name) == 0) return;

        // --- Cluster-level correlations ---
        std::unordered_map<uint64_t, std::unordered_map<std::string, std::vector<std::shared_ptr<Cluster>>>> events;
        for (auto const& [det, clusters] : batch.clusters) {
            for (auto const& cl : clusters) events[cl->eventID()][det].push_back(cl);
        }

        std::vector<decltype(events)::value_type*> event_ptrs;
        event_ptrs.reserve(events.size());
        for (auto& kv : events) event_ptrs.push_back(&kv);

        std::vector<std::unordered_map<std::string, DetAccum>> thread_accums(n_threads);
        size_t chunk = (event_ptrs.size() + n_threads - 1) / std::max<size_t>(n_threads, 1);
        std::vector<std::future<void>> futures;

        for (size_t t = 0; t < n_threads; ++t) {
            size_t start = t * chunk, end = std::min(start + chunk, event_ptrs.size());
            if (start >= end) continue;
            futures.push_back(thread_pool->submit([&, start, end, t]() {
                auto& local = thread_accums[t];
                for (size_t k = start; k < end; ++k) {
                    auto& dets = event_ptrs[k]->second;
                    if (dets.count(ref_name) == 0) continue;
                    auto const& ref_cls = dets.at(ref_name);

                    for (auto const& [name, clusters] : dets) {
                        std::string dir = getName() + "/" + name;
                        if (!pm.hasPlot1D(dir, "correlationX")) continue;

                        auto& det = geo.getDetector(name);
                        auto& ref_det = geo.getDetector(ref_name);
                        auto& acc = local[name];
                        if (acc.corrX2Dlocal.size() == 0) {
                            acc.corrX2Dlocal = Eigen::MatrixXd::Zero(det.n_pixels_x, ref_det.n_pixels_x);
                            acc.corrY2Dlocal = Eigen::MatrixXd::Zero(det.n_pixels_y, ref_det.n_pixels_y);
                            acc.corrX2D = Eigen::MatrixXd::Zero(100, 100);
                            acc.corrY2D = Eigen::MatrixXd::Zero(100, 100);
                            acc.corrXY2D = Eigen::MatrixXd::Zero(100, 100);
                            acc.corrYX2D = Eigen::MatrixXd::Zero(100, 100);
                        }

                        for (auto const& c : clusters) {
                            for (auto const& r : ref_cls) {
                                double time_diff = r->timestamp() - c->timestamp();

                                if (std::abs(time_diff) < time_cut_abs_ || !do_time_cut_) {
                                    double dx = r->global().x() - c->global().x();
                                    double dy = r->global().y() - c->global().y();

                                    acc.corrX.push_back(dx);
                                    acc.corrY.push_back(dy);
                                    // Matches corry's own (seemingly crossed
                                    // but deliberate) naming exactly: "XY"
                                    // pairs the REFERENCE's Y with THIS
                                    // detector's X; "YX" pairs the
                                    // reference's X with this detector's Y.
                                    acc.corrXY.push_back(r->global().y() - c->global().x());
                                    acc.corrYX.push_back(r->global().x() - c->global().y());

                                    fill2DLocal(acc.corrX2Dlocal, acc.momX2Dlocal, c->column(), r->column(), -0.5, det.n_pixels_x - 0.5, -0.5, ref_det.n_pixels_x - 0.5);
                                    fill2DLocal(acc.corrY2Dlocal, acc.momY2Dlocal, c->row(), r->row(), -0.5, det.n_pixels_y - 0.5, -0.5, ref_det.n_pixels_y - 0.5);

                                    fill2DLocal(acc.corrX2D, acc.momX2D, c->global().x(), r->global().x(), -10.1, 9.9, -10.1, 9.9);
                                    fill2DLocal(acc.corrY2D, acc.momY2D, c->global().y(), r->global().y(), -10.1, 9.9, -10.1, 9.9);
                                    fill2DLocal(acc.corrXY2D, acc.momXY2D, r->global().y(), c->global().x(), -10.1, 9.9, -10.1, 9.9);
                                    fill2DLocal(acc.corrYX2D, acc.momYX2D, r->global().x(), c->global().y(), -10.1, 9.9, -10.1, 9.9);
                                }

                                acc.corrTime.push_back(time_diff);
                                // 40 MHz clock ticks (25 ns each), matching corry.
                                acc.corrTimeInt.push_back(static_cast<double>(static_cast<long long>(time_diff / 25)));
                            }
                        }
                    }
                }
            }));
        }
        for (auto& f : futures) f.get();

        std::unordered_map<std::string, DetAccum> merged;
        for (auto& local : thread_accums) {
            for (auto& [name, acc] : local) {
                auto& m = merged[name];
                if (m.corrX2Dlocal.size() == 0) {
                    m.corrX2Dlocal = acc.corrX2Dlocal; m.corrY2Dlocal = acc.corrY2Dlocal;
                    m.corrX2D = acc.corrX2D; m.corrY2D = acc.corrY2D; m.corrXY2D = acc.corrXY2D; m.corrYX2D = acc.corrYX2D;
                } else {
                    m.corrX2Dlocal += acc.corrX2Dlocal; m.corrY2Dlocal += acc.corrY2Dlocal;
                    m.corrX2D += acc.corrX2D; m.corrY2D += acc.corrY2D; m.corrXY2D += acc.corrXY2D; m.corrYX2D += acc.corrYX2D;
                }
                m.momX2Dlocal += acc.momX2Dlocal; m.momY2Dlocal += acc.momY2Dlocal;
                m.momX2D += acc.momX2D; m.momY2D += acc.momY2D; m.momXY2D += acc.momXY2D; m.momYX2D += acc.momYX2D;
                m.corrX.insert(m.corrX.end(), acc.corrX.begin(), acc.corrX.end());
                m.corrY.insert(m.corrY.end(), acc.corrY.begin(), acc.corrY.end());
                m.corrXY.insert(m.corrXY.end(), acc.corrXY.begin(), acc.corrXY.end());
                m.corrYX.insert(m.corrYX.end(), acc.corrYX.begin(), acc.corrYX.end());
                m.corrTime.insert(m.corrTime.end(), acc.corrTime.begin(), acc.corrTime.end());
                m.corrTimeInt.insert(m.corrTimeInt.end(), acc.corrTimeInt.begin(), acc.corrTimeInt.end());
            }
        }
        for (auto& [name, acc] : merged) {
            std::string dir = getName() + "/" + name;
            pm.getPlot2D(dir, "correlationX_2Dlocal").fillBatch(acc.corrX2Dlocal, acc.momX2Dlocal.sw, acc.momX2Dlocal.sx,
                                acc.momX2Dlocal.sy, acc.momX2Dlocal.sx2, acc.momX2Dlocal.sy2, acc.momX2Dlocal.sxy);
            pm.getPlot2D(dir, "correlationY_2Dlocal").fillBatch(acc.corrY2Dlocal, acc.momY2Dlocal.sw, acc.momY2Dlocal.sx,
                                acc.momY2Dlocal.sy, acc.momY2Dlocal.sx2, acc.momY2Dlocal.sy2, acc.momY2Dlocal.sxy);
            pm.getPlot2D(dir, "correlationX_2D").fillBatch(acc.corrX2D, acc.momX2D.sw, acc.momX2D.sx,
                                acc.momX2D.sy, acc.momX2D.sx2, acc.momX2D.sy2, acc.momX2D.sxy);
            pm.getPlot2D(dir, "correlationY_2D").fillBatch(acc.corrY2D, acc.momY2D.sw, acc.momY2D.sx,
                                acc.momY2D.sy, acc.momY2D.sx2, acc.momY2D.sy2, acc.momY2D.sxy);
            pm.getPlot2D(dir, "correlationXY_2D").fillBatch(acc.corrXY2D, acc.momXY2D.sw, acc.momXY2D.sx,
                                acc.momXY2D.sy, acc.momXY2D.sx2, acc.momXY2D.sy2, acc.momXY2D.sxy);
            pm.getPlot2D(dir, "correlationYX_2D").fillBatch(acc.corrYX2D, acc.momYX2D.sw, acc.momYX2D.sx,
                                acc.momYX2D.sy, acc.momYX2D.sx2, acc.momYX2D.sy2, acc.momYX2D.sxy);
            pm.getPlot1D(dir, "correlationX").fillBatch(acc.corrX);
            pm.getPlot1D(dir, "correlationY").fillBatch(acc.corrY);
            pm.getPlot1D(dir, "correlationXY").fillBatch(acc.corrXY);
            pm.getPlot1D(dir, "correlationYX").fillBatch(acc.corrYX);
            pm.getPlot1D(dir, "correlationTime").fillBatch(acc.corrTime);
            pm.getPlot1D(dir, "correlationTimeInt").fillBatch(acc.corrTimeInt);
        }
    }

    REGISTER_MODULE(Correlations)
}
