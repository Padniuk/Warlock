#include "WaveformSelector.hpp"
#include "core/ModuleFactory.hpp"
#include "core/utils/PlotManager.hpp"
#include "core/utils/Logger.hpp"
#include "core/detector/GeometryManager.hpp"
#include <algorithm>

namespace framework {
    WaveformSelector::WaveformSelector(Configuration cfg, Configuration g_cfg, ThreadPool* pool)
        : Module(std::move(cfg), std::move(g_cfg), pool) {
        detector_name_ = config.get<std::string>("name", "");

        min_snr_ = config.get<double>("min_snr", 0.0);
        max_snr_ = config.get<double>("max_snr", 1e30);
        min_amplitude_ = config.get<double>("min_amplitude", 0.0);
        max_amplitude_ = config.get<double>("max_amplitude", 1e30);
        min_noise_ = config.get<double>("min_noise", 0.0);
        max_noise_ = config.get<double>("max_noise", 1e30);
        min_rise_time_ = config.get<double>("min_rise_time", 0.0);
        max_rise_time_ = config.get<double>("max_rise_time", 1e30);
        min_start_time_ = config.get<double>("min_start_time", 0.0);
        max_start_time_ = config.get<double>("max_start_time", 1e30);
        min_integral_ = config.get<double>("min_integral", 0.0);
        max_integral_ = config.get<double>("max_integral", 1e30);
        min_time_over_threshold_ = config.get<double>("min_time_over_threshold", 0.0);
        max_time_over_threshold_ = config.get<double>("max_time_over_threshold", 1e30);
        min_time_over_noise_ = config.get<double>("min_time_over_noise", 0.0);
        max_time_over_noise_ = config.get<double>("max_time_over_noise", 1e30);

        cfd_ = config.get<double>("cfd", 0.5);

        dut_names_ = config.getArray<std::string>("dut_names");
    }

    void WaveformSelector::registerChannelPlots(const std::string& dir) {
        auto& pm = PlotManager::getInstance();
        pm.registerPlot1D(dir, "snr", 500, 0, 200);
        pm.registerPlot1D(dir, "amplitude", 500, 0, 1.2);
        pm.registerPlot1D(dir, "rise_time", 200, 0, 5.0);
        pm.registerPlot1D(dir, "baseline", 500, -0.1, 0.1);
        pm.registerPlot1D(dir, "noise", 500, 0, 0.05);
        pm.registerPlot1D(dir, "start_time", 500, 0, 500);
        pm.registerPlot1D(dir, "integral", 500, 0, 50.0);
        pm.registerPlot1D(dir, "time_over_threshold", 500, 0, 50);
        pm.registerPlot1D(dir, "time_over_noise", 500, 0, 200);
    }

    void WaveformSelector::fillChannelPlots(const std::string& dir, double snr, double amplitude, double noise,
                                             double baseline, double rise_time, double start_time, double integral,
                                             double time_over_threshold, double time_over_noise) {
        auto& pm = PlotManager::getInstance();
        pm.fill1D(dir, "snr", snr);
        pm.fill1D(dir, "amplitude", amplitude);
        pm.fill1D(dir, "rise_time", rise_time);
        pm.fill1D(dir, "baseline", baseline);
        pm.fill1D(dir, "noise", noise);
        pm.fill1D(dir, "start_time", start_time);
        pm.fill1D(dir, "integral", integral);
        pm.fill1D(dir, "time_over_threshold", time_over_threshold);
        pm.fill1D(dir, "time_over_noise", time_over_noise);
    }

    void WaveformSelector::initialize() {
        WR_LOG(STATUS, "Initializing Waveform Selector...");
        auto& geo = GeometryManager::getInstance();
        if (dut_names_.empty()) {
            if (!detector_name_.empty()) {
                // Named instance (one WaveformSelector per detector, via
                // repeated [[WaveformSelector]] blocks): only this
                // detector's histograms are relevant.
                dut_names_.push_back(detector_name_);
            } else {
                for (const auto& name : geo.getDetectorNames()) {
                    if (geo.getDetector(name).type == "caendt5742") dut_names_.push_back(name);
                }
            }
        }

        // Same output shape as WaveformProcessingCAEN (total + per-detector
        // total/channel breakdown), but without a CFD level - the loaded
        // waveforms only carry a single, already-fixed start_time (whichever
        // CFD fraction was used when they were originally processed and
        // saved), so there's no CFD axis left to break down at this stage.
        registerChannelPlots(getName() + "/total");
        for (const auto& det_name : dut_names_) {
            if (!geo.hasDetector(det_name)) continue;
            auto& det = geo.getDetector(det_name);

            std::string det_dir = getName() + "/" + det_name;
            registerChannelPlots(det_dir + "/total");

            for (int row = 0; row < det.n_pixels_y; ++row) {
                for (int col = 0; col < det.n_pixels_x; ++col) {
                    int channel = col + det.n_pixels_x * row;
                    registerChannelPlots(det_dir + "/channel_" + std::to_string(channel));
                }
            }
        }
    }

    void WaveformSelector::makePixel(DataBatch& batch, const std::string& det_name, const Waveform& wf) const {
        auto& geo = GeometryManager::getInstance();
        if (!geo.hasDetector(det_name)) return;

        // charge = the computed signal integral; timestamp = event baseline
        // time plus the signal's own start_time offset (at this instance's
        // configured CFD fraction) - same convention as Warlock's original
        // (pre-rewrite) WaveformSelector prototype. ClusteringSpatial merges
        // this with any adjacent-channel pixel from the same event.
        auto pixel = std::make_shared<Pixel>(det_name, wf.column(), wf.row(), wf.integral,
                                              wf.timestamp() + wf.startTimeAtCfd(cfd_), wf.eventID(), wf.channel());
        batch.pixels[det_name].push_back(pixel);
    }

    bool WaveformSelector::passesCuts(const Waveform& wf) const {
        double start_time = wf.startTimeAtCfd(cfd_);
        return wf.snr >= min_snr_ && wf.snr <= max_snr_ &&
               wf.amplitude >= min_amplitude_ && wf.amplitude <= max_amplitude_ &&
               wf.noise >= min_noise_ && wf.noise <= max_noise_ &&
               wf.rise_time >= min_rise_time_ && wf.rise_time <= max_rise_time_ &&
               start_time >= min_start_time_ && start_time <= max_start_time_ &&
               wf.integral >= min_integral_ && wf.integral <= max_integral_ &&
               wf.time_over_threshold >= min_time_over_threshold_ && wf.time_over_threshold <= max_time_over_threshold_ &&
               wf.time_over_noise >= min_time_over_noise_ && wf.time_over_noise <= max_time_over_noise_;
    }

    void WaveformSelector::run(DataBatch& batch) {
        auto& geo = GeometryManager::getInstance();
        for (auto& entry : batch.waveforms) {
            const std::string& det_name = entry.first;

            // A named instance only touches its own detector - every other
            // key in batch.waveforms is left completely alone, so several
            // [[WaveformSelector]] instances (one per detector, each with
            // its own cuts) can be chained without stepping on each other.
            if (!detector_name_.empty() && det_name != detector_name_) continue;

            auto& waveforms = entry.second;

            bool is_dut = std::find(dut_names_.begin(), dut_names_.end(), det_name) != dut_names_.end();
            int n_pixels_x = geo.hasDetector(det_name) ? geo.getDetector(det_name).n_pixels_x : 0;

            std::vector<std::shared_ptr<Waveform>> accepted;
            accepted.reserve(waveforms.size());

            for (auto& wf : waveforms) {
                if (!wf->is_processed) continue;
                total_in_++;

                if (!passesCuts(*wf)) continue;

                double start_time = wf->startTimeAtCfd(cfd_);
                fillChannelPlots(getName() + "/total", wf->snr, wf->amplitude, wf->noise,
                                  wf->baseline, wf->rise_time, start_time, wf->integral,
                                  wf->time_over_threshold, wf->time_over_noise);

                if (is_dut) {
                    std::string det_dir = getName() + "/" + det_name;
                    int channel = n_pixels_x > 0 ? (wf->column() + n_pixels_x * wf->row()) : 0;
                    fillChannelPlots(det_dir + "/total", wf->snr, wf->amplitude, wf->noise,
                                      wf->baseline, wf->rise_time, start_time, wf->integral,
                                      wf->time_over_threshold, wf->time_over_noise);
                    fillChannelPlots(det_dir + "/channel_" + std::to_string(channel), wf->snr, wf->amplitude,
                                      wf->noise, wf->baseline, wf->rise_time, start_time, wf->integral,
                                      wf->time_over_threshold, wf->time_over_noise);
                }

                total_out_++;
                accepted.push_back(wf);
                makePixel(batch, det_name, *wf);
            }

            // Replace in place - downstream modules (DUT association,
            // analysis efficiency, analysis timing) read batch.waveforms
            // and should only ever see waveforms that passed selection.
            waveforms = std::move(accepted);
        }
    }

    void WaveformSelector::finalize() {
        WR_LOG(STATUS, "Waveform Selector: " + std::to_string(total_out_) + " / " + std::to_string(total_in_) + " waveforms passed.");
    }
    REGISTER_MODULE(WaveformSelector)
}
