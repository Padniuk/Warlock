#include "ClusteringSpatial.hpp"
#include "core/ModuleFactory.hpp"
#include "core/detector/GeometryManager.hpp"
#include "core/utils/PlotManager.hpp"
#include "core/utils/Logger.hpp"
#include <set>
#include <future>
#include <map>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <sstream>

namespace framework {

    namespace {
        // CPU-time summed across every chunk task in every batch - compared
        // against the module's own wall-clock time (FrameworkManager's
        // "Wall-clock timing" summary) to gauge how much the per-detector
        // event-range chunking (see run() below) actually parallelizes vs.
        // losing time back to per-detector histogram lock contention
        // (multiple chunks of the SAME detector fill the SAME "clusterSize"
        // etc. histograms concurrently). Plain double: summed from
        // per-chunk local values after every future in a batch is joined,
        // single-threaded at that point.
        double stat_time_total_s = 0.0;
    }

    ClusteringSpatial::ClusteringSpatial(Configuration cfg, Configuration g_cfg, ThreadPool* pool)
        : Module(std::move(cfg), std::move(g_cfg), pool) {
        use_weighting_ = config.get<bool>("charge_weighting", true);
        discard_events_missing_on_ = config.getArray<std::string>("discard_events_missing_on");
        // Accepts either detector names or detector types, as a single
        // bare value or a list - same exclude-list convention as
        // Tracking4D::exclude_detectors_/AlignmentMillepede::
        // exclude_detectors_/EventLoaderHDF5::exclude_detectors_.
        exclude_detectors_ = config.getArray<std::string>("exclude_detectors");
    }

    bool ClusteringSpatial::isExcluded(const DetectorGeo& det) const {
        return std::find(exclude_detectors_.begin(), exclude_detectors_.end(), det.name) != exclude_detectors_.end() ||
               std::find(exclude_detectors_.begin(), exclude_detectors_.end(), det.type) != exclude_detectors_.end();
    }

    void ClusteringSpatial::initialize() {
        WR_LOG(STATUS, "Initializing ClusteringSpatial...");
        auto& geo = GeometryManager::getInstance();
        
        for (const auto& name : geo.getDetectorNames()) {
            auto& det = geo.getDetector(name);
            if (!det.isPixelDetector()) continue;

            // Masked pixels are now loaded once by GeometryManager itself
            // (shared across every module that needs them) instead of each
            // module parsing its own copy - see det.masked_pixels below.
            if (!det.masked_pixels.empty()) {
                WR_LOG(INFO, "Loaded " + std::to_string(det.masked_pixels.size()) + " masked pixels for " + name);
            }

            std::string dir = getName() + "/" + name;
            auto& pm = PlotManager::getInstance();

            // Binning matched exactly to corryvreckan's ClusteringSpatial.
            pm.registerPlot1D(dir, "clusterSize", 100, -0.5, 99.5);
            pm.registerPlot1D(dir, "clusterMultiplicity", 50, -0.5, 49.5);
            pm.registerPlot1D(dir, "clusterSeedCharge", 256, -0.5, 255.5);
            pm.registerPlot1D(dir, "clusterCharge", 5000, -0.5, 49999.5);
            pm.registerPlot1D(dir, "clusterWidthRow", 25, -0.5, 24.5);
            pm.registerPlot1D(dir, "clusterWidthColumn", 100, -0.5, 99.5);
            pm.registerPlot1D(dir, "clusterTimes", 10000, 0, 100);
            // corry's uncertainty plots are in um, ranged to the pitch.
            pm.registerPlot1D(dir, "clusterUncertaintyX", 100, 0, det.pitch_x * 1000.0);
            pm.registerPlot1D(dir, "clusterUncertaintyY", 100, 0, det.pitch_y * 1000.0);

            pm.registerPlot2D(dir, "clusterPositionLocal", det.n_pixels_x, -0.5, det.n_pixels_x - 0.5, det.n_pixels_y, -0.5, det.n_pixels_y - 0.5);
            // corry centers this on the detector's global position and spans
            // its physical size, rather than a fixed generic window.
            double half_size_x = det.sizeX() / 2.0;
            double half_size_y = det.sizeY() / 2.0;
            pm.registerPlot2D(dir, "clusterPositionGlobal", 400, det.position.x() - half_size_x, det.position.x() + half_size_x,
                                                             400, det.position.y() - half_size_y, det.position.y() + half_size_y);
        }
    }

    void ClusteringSpatial::run(DataBatch& batch) {
        if (batch.pixels.empty()) return;
        std::vector<std::future<void>> futures;

        std::set<uint64_t> batch_events;
        if (discard_events_missing_on_.empty()) {
            for (auto const& entry : batch.pixels) {
                for (auto const& p : entry.second) batch_events.insert(p->eventID());
            }
        } else {
            bool first = true;
            for (const auto& req_det : discard_events_missing_on_) {
                std::set<uint64_t> det_events;
                if (batch.pixels.count(req_det)) {
                    for (auto const& p : batch.pixels.at(req_det)) det_events.insert(p->eventID());
                }
                if (first) {
                    batch_events = det_events;
                    first = false;
                } else {
                    std::set<uint64_t> intersect;
                    std::set_intersection(batch_events.begin(), batch_events.end(),
                                          det_events.begin(), det_events.end(),
                                          std::inserter(intersect, intersect.begin()));
                    batch_events = intersect;
                }
            }
        }

        auto& geo_filter = GeometryManager::getInstance();
        for (auto& entry : batch.pixels) {
            const std::string det_name = entry.first;
            if (geo_filter.hasDetector(det_name) && isExcluded(geo_filter.getDetector(det_name))) continue;
            auto& pixels = entry.second;

            // This detector is about to be (re)clustered from its own
            // pixels this batch - discard whatever clusters it already has
            // (e.g. loaded from an earlier stage's saved output, before a
            // fresh clustering pass over the same detector's raw waveforms
            // in this config) rather than letting the two silently coexist.
            // Stale clusters left in place would give downstream modules
            // like DUTAssociation two independently-computed candidates per
            // track and split timing/efficiency plots across both. The
            // per-chunk `insert()` below (merging this detector's own
            // concurrently-processing thread chunks together) is correct
            // and left alone - this clear() only removes clusters from
            // before this detector's clustering pass even started.
            batch.clusters[det_name].clear();

            // Built once per detector, read-only from here on - shared via
            // shared_ptr so every sub-chunk task below (see next comment)
            // can safely hold a reference to it after this loop iteration's
            // stack frame is gone.
            auto event_map = std::make_shared<std::map<uint64_t, std::vector<std::shared_ptr<Pixel>>>>();
            for (auto& p : pixels) (*event_map)[p->eventID()].push_back(p);

            // A single task per detector caps usable parallelism at the
            // pixel-detector count regardless of the configured thread
            // count, so each detector's own event range is further split
            // into n_threads sub-chunks (the same pattern Tracking4D/
            // WaveformProcessingCAEN use) to spread every detector's work
            // across all available threads.
            //
            // events_vec is a shared_ptr, not a plain local vector captured
            // by reference: it's loop-local (rebuilt fresh each detector
            // iteration), unlike `pixels` above which refers into
            // batch.pixels' own storage and outlives the whole function - a
            // queued task for THIS detector can still be waiting to run (or
            // running) once the loop moves to the next detector and this
            // variable's stack slot is reused, so capturing it by reference
            // is a dangling-reference bug that silently corrupts which
            // events get clustered.
            auto events_vec = std::make_shared<std::vector<uint64_t>>(batch_events.begin(), batch_events.end());
            size_t n_threads = thread_pool->getThreadCount();
            size_t chunk_size = (events_vec->size() + n_threads - 1) / n_threads;
            if (chunk_size == 0) chunk_size = 1;

            for (size_t chunk_start = 0; chunk_start < events_vec->size(); chunk_start += chunk_size) {
            size_t chunk_end = std::min(chunk_start + chunk_size, events_vec->size());

            futures.push_back(thread_pool->submit([this, det_name, chunk_start, chunk_end, events_vec, event_map, &batch]() {
                auto t_task_start = std::chrono::steady_clock::now();
                auto& geo = GeometryManager::getInstance();
                auto& det = geo.getDetector(det_name);
                auto& pm = PlotManager::getInstance();
                std::string dir = getName() + "/" + det_name;

                bool has_mask = !det.masked_pixels.empty();

                std::vector<std::shared_ptr<Cluster>> detector_batch_clusters;

                for (size_t ev_idx = chunk_start; ev_idx < chunk_end; ++ev_idx) {
                    uint64_t ev_id = (*events_vec)[ev_idx];
                    size_t clusters_in_event = 0;

                    auto it_ev = event_map->find(ev_id);
                    if (it_ev != event_map->end()) {
                        auto& ev_pixels = it_ev->second;
                        std::map<std::pair<int, int>, std::shared_ptr<Pixel>> hitmap;
                        
                        for (auto& p : ev_pixels) {
                            if (has_mask && det.isMasked(p->column(), p->row())) continue;
                            hitmap[{p->column(), p->row()}] = p;
                        }

                        std::set<std::shared_ptr<Pixel>> used_pixels;
                        
                        for (auto& [coords, p] : hitmap) {
                            if (used_pixels.count(p)) continue;
                            
                            auto cluster = std::make_shared<Cluster>();
                            cluster->setDetector(&det);
                            cluster->setEventID(ev_id);

                            std::vector<std::shared_ptr<Pixel>> queue = {p};
                            used_pixels.insert(p);
                            size_t head = 0;
                            
                            while(head < queue.size()){
                                auto curr = queue[head++];
                                cluster->addPixel(curr);
                                
                                for(int dx = -1; dx <= 1; ++dx){
                                    for(int dy = -1; dy <= 1; ++dy){
                                        if(dx == 0 && dy == 0) continue;
                                        auto it = hitmap.find({curr->column() + dx, curr->row() + dy});
                                        if(it != hitmap.end() && !used_pixels.count(it->second)){
                                            used_pixels.insert(it->second);
                                            queue.push_back(it->second);
                                        }
                                    }
                                }
                            }

                            double sum_col = 0, sum_row = 0, total_q = 0, sum_ts = 0;
                            double seed_q = 0;
                            int min_col = 1e9, max_col = -1e9, min_row = 1e9, max_row = -1e9;

                            for (auto& cp : cluster->pixels()) {
                                double weight = use_weighting_ ? cp->charge() : 1.0;
                                sum_col += cp->column() * weight;
                                sum_row += cp->row() * weight;
                                sum_ts += cp->timestamp() * weight;
                                total_q += weight;

                                if (cp->charge() > seed_q) seed_q = cp->charge();
                                if (cp->column() < min_col) min_col = cp->column();
                                if (cp->column() > max_col) max_col = cp->column();
                                if (cp->row() < min_row) min_row = cp->row();
                                if (cp->row() > max_row) max_row = cp->row();
                            }
                            
                            double cog_col = sum_col / total_q;
                            double cog_row = sum_row / total_q;

                            cluster->setTimestamp(sum_ts / total_q);
                            cluster->setLocalCenter(cog_col, cog_row);

                            double lx = det.getLocalX(cog_col);
                            double ly = det.getLocalY(cog_row);
                            
                            Eigen::AngleAxisd rotX(det.orientation.x(), Eigen::Vector3d::UnitX());
                            Eigen::AngleAxisd rotY(det.orientation.y(), Eigen::Vector3d::UnitY());
                            Eigen::AngleAxisd rotZ(det.orientation.z(), Eigen::Vector3d::UnitZ());
                            
                            Eigen::Matrix3d R = (rotZ * rotY * rotX).matrix();
                            Eigen::Vector3d local_pos(lx, ly, 0.0);
                            Eigen::Vector3d global_pos = R * local_pos + det.position;
                            
                            cluster->setGlobal(global_pos);

                            pm.fill1D(dir, "clusterSize", static_cast<double>(cluster->size()));
                            pm.fill1D(dir, "clusterSeedCharge", seed_q);
                            pm.fill1D(dir, "clusterCharge", total_q);
                            pm.fill1D(dir, "clusterWidthColumn", static_cast<double>(max_col - min_col + 1));
                            pm.fill1D(dir, "clusterWidthRow", static_cast<double>(max_row - min_row + 1));
                            pm.fill1D(dir, "clusterTimes", cluster->timestamp());
                            // corry's uncertainty plots are filled in um, matching the
                            // registered pitch-scaled range above.
                            pm.fill1D(dir, "clusterUncertaintyX", det.spatial_res_x * 1000.0);
                            pm.fill1D(dir, "clusterUncertaintyY", det.spatial_res_y * 1000.0);
                            pm.fill2D(dir, "clusterPositionLocal", cog_col, cog_row);
                            pm.fill2D(dir, "clusterPositionGlobal", global_pos.x(), global_pos.y());

                            detector_batch_clusters.push_back(std::move(cluster));
                            clusters_in_event++;
                        }
                        // Matches corry: ClusteringSpatial::run() returns
                        // early (before ever reaching its own
                        // clusterMultiplicity fill) whenever THIS detector
                        // has zero pixels for the current event, rather
                        // than filling a 0 - batch_events is the union of
                        // every detector's events, so an event can be "in
                        // scope" here (some OTHER detector had a hit) while
                        // this one had none.
                        pm.fill1D(dir, "clusterMultiplicity", static_cast<double>(clusters_in_event));
                    }
                }

                double task_time = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_task_start).count();
                {
                    std::lock_guard<std::mutex> lock(cluster_mtx_);
                    batch.clusters[det_name].insert(batch.clusters[det_name].end(),
                                                   detector_batch_clusters.begin(), detector_batch_clusters.end());
                    stat_time_total_s += task_time;
                }
            }));
            }
        }

        for (auto& f : futures) f.get();
    }

    void ClusteringSpatial::finalize() {
        std::ostringstream stats;
        stats << "ClusteringSpatial CPU-time (seconds, summed across worker threads): total=" << stat_time_total_s;
        WR_LOG(STATUS, stats.str());
    }

    REGISTER_MODULE(ClusteringSpatial)
}