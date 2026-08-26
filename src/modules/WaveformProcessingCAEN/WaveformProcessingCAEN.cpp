/** @file WaveformProcessingCAEN.cpp */
#include "WaveformProcessingCAEN.hpp"
#include "core/ModuleFactory.hpp"
#include "core/utils/PlotManager.hpp"
#include "core/utils/Logger.hpp"
#include "core/detector/GeometryManager.hpp"
#include <Eigen/Dense>
#include <future>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <mutex>
#include <sstream>
#include <iomanip>
#include <vector>
#include <chrono>

namespace framework {

    namespace {
        /// CAEN DT5742 ADC-code -> Volts calibration constants, mirrored
        /// from RawConverter.cpp (kCaenAdcMid/kCaenAdcLsbDen/kCaenVpp) and,
        /// beneath that, eudaq's own CAENDT5742RawEvent2StdEventConverter.cc
        /// ADC12ToVolts(). "/hits/samples" now stores the DRS4-corrected ADC
        /// code (int16), not calibrated Volts - minionRawConverter defers
        /// this last, cheap, stateless step to here deliberately, since it
        /// only needs paying for on whatever fraction of conversions actually
        /// load waveforms, unlike DRS4 correction itself (expensive,
        /// stateful, correction-table-dependent - has to stay in the
        /// converter). Keep these three constants in sync with
        /// RawConverter.cpp's copy if either ever changes.
        constexpr double kCaenAdcMid = 2048.0;
        constexpr double kCaenAdcLsbDen = 4096.0;
        constexpr double kCaenVpp = 1.0;

        /// Cumulative CPU-time (sum across all worker threads, not
        /// wall-clock) spent inside plot-fill calls vs everything else per
        /// waveform - unlike Tracking4D, this module already fills plots
        /// from inside its parallel workers, so there's no serial phase to
        /// find; the open question is whether the SHARED "raw_*" histograms
        /// (registered once under getName(), hit by every waveform from
        /// every thread) cause lock contention expensive enough to matter.
        /// Plain double, summed from per-chunk local accumulators after
        /// every future in a batch has been joined (single-threaded at that
        /// point) - never touched concurrently itself, so no atomic needed.
        double stat_time_fill_s = 0.0;
        double stat_time_math_s = 0.0;

        /// RAII wall-clock accumulator: added to `accum` on scope exit
        /// regardless of which `continue` path leaves the scope, so a timer
        /// placed at the top of a loop body with early-exit branches still
        /// captures every iteration's full duration correctly.
        struct ScopeTimer {
            std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
            double& accum;
            explicit ScopeTimer(double& a) : accum(a) {}
            ~ScopeTimer() { accum += std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count(); }
        };
        /// Ported from corry's WfSignal.cpp `get_median()` - median via
        /// nth_element (partial sort), including its exact even/odd-count
        /// averaging, so this matches corry's own median-based signal-start
        /// detection bit-for-bit rather than just "a" median implementation.
        /// Takes the vector by reference and reorders it in place (nth_element
        /// always does this regardless) rather than by value: both call sites
        /// below only ever use the vector to compute this one median and
        /// never need its original order again afterward.
        double medianOf(std::vector<double>& values) {
            if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
            size_t n = values.size();
            size_t mid = n / 2;
            std::nth_element(values.begin(), values.begin() + static_cast<long>(mid), values.end());
            double median = values[mid];
            if (n % 2 == 0) {
                double first = median;
                std::nth_element(values.begin(), values.begin() + static_cast<long>(mid) - 1, values.end());
                median = (first + values[mid - 1]) / 2.0;
            }
            return median;
        }

        /// Linear interpolation for the time at which the rising edge
        /// crosses `target_v`, given `pts` (voltage, time) SORTED ascending
        /// by voltage - the free-standing half of the mechanism ported from
        /// corry's WfPositiveSignal::findTimeAtRisingEdge() (WfSignal.cpp),
        /// which builds a ROOT TGraph(voltage, time) over every sample from
        /// where the pulse leaves its 99%-of-amplitude peak down to the
        /// threshold being searched for, then calls TGraph::Eval(target_v) -
        /// plain (unsplined) TGraph::Eval is itself just a bracket-and-
        /// linearly-interpolate lookup, so this reproduces that same
        /// result without linking ROOT into Warlock's core (which never
        /// otherwise depends on it - see the top-level README). Clamps to
        /// the nearest edge if `target_v` falls outside `pts`' own range,
        /// which by construction of the caller's search shouldn't happen
        /// in practice - only a defined fallback for the boundary case.
        double interpolateTimeAtVoltage(const std::vector<std::pair<double, double>>& pts, double target_v) {
            if (pts.empty()) return 0.0;
            if (pts.size() == 1 || target_v <= pts.front().first) return pts.front().second;
            if (target_v >= pts.back().first) return pts.back().second;
            auto it = std::lower_bound(pts.begin(), pts.end(), target_v,
                [](const std::pair<double, double>& p, double v) { return p.first < v; });
            const auto& hi = *it;
            const auto& lo = *(it - 1);
            double dv = hi.first - lo.first;
            if (dv == 0.0) return lo.second;
            double frac = (target_v - lo.first) / dv;
            return lo.second + frac * (hi.second - lo.second);
        }
    }
    WaveformProcessingCAEN::WaveformProcessingCAEN(Configuration cfg, Configuration g_cfg, ThreadPool* pool)
        : Module(std::move(cfg), std::move(g_cfg), pool) {
        polarity_ = config.get<int>("signal_polarity", -1);

        // `cfd` may be a single scalar or an array in the TOML - normalize to
        // a vector either way, falling back to corry's default of 0.7.
        cfd_values_ = config.getArray<double>("cfd");
        if (cfd_values_.empty()) cfd_values_ = {config.get<double>("cfd", 0.7)};

        sampling_interval_ = config.get<double>("sampling_interval", 0.2);
        noise_factor_ = config.get<double>("noise_factor", 1.0); // corry's own default is 1, not 3
        integral_time_ = config.get<double>("integral_time", 2.8);

        // Matches corry's own config defaults exactly (WaveformProcessingCAEN.cpp
        // constructor) - permissive by default, but still strict-positive on
        // every quantity (min_* = 0.0), so a waveform with no valid signal
        // (snr/start_time/rise_time/integral landing at exactly 0) is
        // rejected even with no explicit cuts configured, same as corry.
        min_snr_ = config.get<double>("min_snr", 0.0);
        max_snr_ = config.get<double>("max_snr", 1e30);
        min_start_time_ = config.get<double>("min_start_time", 0.0);
        max_start_time_ = config.get<double>("max_start_time", 1e30);
        min_rise_time_ = config.get<double>("min_rise_time", 0.0);
        max_rise_time_ = config.get<double>("max_rise_time", 1e30);
        min_integral_ = config.get<double>("min_integral", 0.0);
        max_integral_ = config.get<double>("max_integral", 1e30);

        // Detectors to break down per-channel. Empty means "every
        // caendt5742 detector in the geometry file" (resolved in
        // initialize(), once the geometry is guaranteed to be loaded).
        dut_names_ = config.getArray<std::string>("dut_names");

        storage_file_ = config.get<std::string>("storage_file", "");
        save_enabled_ = !storage_file_.empty();
    }

    std::string WaveformProcessingCAEN::cfdDirName(double cfd) {
        std::ostringstream oss;
        oss << "cfd_" << std::fixed << std::setprecision(2) << cfd;
        return oss.str();
    }

    void WaveformProcessingCAEN::registerChannelPlots(const std::string& dir) {
        auto& pm = PlotManager::getInstance();
        // Fixed axes here (unlike corry's own SetCanExtend()-based per-run
        // auto-fit) risk silently truncating real data if a dataset's
        // range exceeds the constant - SNR can reach into the thousands
        // (heavy-tailed) and rise_time can reach several times its typical
        // sub-ns value, so both are widened with real headroom over
        // observed maxima.
        pm.registerPlot1D(dir, "snr", 400, 0, 400);
        pm.registerPlot1D(dir, "amplitude", 500, 0, 1.2);
        pm.registerPlot1D(dir, "rise_time", 500, 0, 50);
        // Real CAEN DT5742 baseline levels sit around 0.4-0.42 in RAW
        // (unflipped) ADC terms - but unlike `amplitude` (clamped to >= 0
        // below), `baseline` is a plain mean of `s`, the raw trace already
        // multiplied by `polarity_` (see run(), `s = raw * polarity_`), so
        // with the default signal_polarity=-1 it comes out around -0.4,
        // not +0.4. Symmetric range so it holds the real value under
        // either polarity sign, since that's a per-run config choice.
        pm.registerPlot1D(dir, "baseline", 500, -1.2, 1.2);
        pm.registerPlot1D(dir, "noise", 500, 0, 0.05);
        pm.registerPlot1D(dir, "start_time", 500, 0, 500);
        // Named "integral" to match the actual quantity/HDF5 dataset name
        // (2026-08-23) - this folder used to be called "charge", which never
        // corresponded to any real per-waveform field (WaveformProcessingCAEN
        // only ever computed/saved "integral"; a downstream analysis script
        // reading this folder as "charge" was silently reading nothing).
        pm.registerPlot1D(dir, "integral", 500, 0, 50.0);
        pm.registerPlot1D(dir, "time_over_threshold", 500, 0, 50);
        pm.registerPlot1D(dir, "time_over_noise", 500, 0, 200);
    }

    WaveformProcessingCAEN::ChannelPlotRefs WaveformProcessingCAEN::fetchChannelPlots(const std::string& dir) {
        auto& pm = PlotManager::getInstance();
        ChannelPlotRefs refs;
        refs.snr = &pm.getPlot1D(dir, "snr");
        refs.amplitude = &pm.getPlot1D(dir, "amplitude");
        refs.noise = &pm.getPlot1D(dir, "noise");
        refs.baseline = &pm.getPlot1D(dir, "baseline");
        refs.rise_time = &pm.getPlot1D(dir, "rise_time");
        refs.start_time = &pm.getPlot1D(dir, "start_time");
        refs.integral = &pm.getPlot1D(dir, "integral");
        refs.time_over_threshold = &pm.getPlot1D(dir, "time_over_threshold");
        refs.time_over_noise = &pm.getPlot1D(dir, "time_over_noise");
        return refs;
    }

    void WaveformProcessingCAEN::fillChannelPlots(const ChannelPlotRefs& refs, double snr, double amplitude, double noise,
                                               double baseline, double rise_time, double start_time, double integral,
                                               double time_over_threshold, double time_over_noise) {
        if (!refs.snr) return; // dir was never registered - matches PlotManager::fill1D()'s no-op contract
        refs.snr->fill(snr);
        refs.amplitude->fill(amplitude);
        refs.rise_time->fill(rise_time);
        refs.baseline->fill(baseline);
        refs.noise->fill(noise);
        refs.start_time->fill(start_time);
        refs.integral->fill(integral);
        refs.time_over_threshold->fill(time_over_threshold);
        refs.time_over_noise->fill(time_over_noise);
    }

    void WaveformProcessingCAEN::titleChannelPlots(const ChannelPlotRefs& refs, const std::string& title) {
        if (!refs.snr) return;
        refs.snr->setTitle(title + " SNR");
        refs.amplitude->setTitle(title + " Amplitude");
        refs.rise_time->setTitle(title + " Rise Time");
        refs.baseline->setTitle(title + " Baseline");
        refs.noise->setTitle(title + " Noise");
        refs.start_time->setTitle(title + " Start Time");
        refs.integral->setTitle(title + " Integral");
        refs.time_over_threshold->setTitle(title + " Time Over Threshold (50%)");
        refs.time_over_noise->setTitle(title + " Time Over Noise");
    }

    void WaveformProcessingCAEN::initialize() {
        WR_LOG(STATUS, "Initializing Waveform Math Engine...");
        auto& pm = PlotManager::getInstance();
        pm.registerPlot1D(getName(), "raw_snr", 500, 0, 5000);
        pm.registerPlot1D(getName(), "raw_amplitude", 500, 0, 1.2);
        pm.registerPlot1D(getName(), "raw_rise_time", 500, 0, 50);
        pm.registerPlot1D(getName(), "raw_baseline", 500, -1.2, 1.2);
        pm.registerPlot1D(getName(), "raw_noise", 500, 0, 0.05);
        pm.registerPlot1D(getName(), "raw_signal_start_time", 500, 0, 500);
        pm.registerPlot1D(getName(), "raw_integral", 500, 0, 50.0);
        pm.registerPlot1D(getName(), "raw_time_over_threshold", 500, 0, 50);
        pm.registerPlot1D(getName(), "raw_time_over_noise", 500, 0, 200);

        auto& geo = GeometryManager::getInstance();
        if (dut_names_.empty()) {
            // Default to the waveform-capable (CAEN digitizer) boards only -
            // NOT every detector in the geometry file. The pixel telescope
            // planes (MIMOSA26, 1152x576 = ~663k pixels each) would each
            // blow up into hundreds of thousands of "channel_<n>" histogram
            // folders if included, exhausting memory almost immediately.
            for (const auto& name : geo.getDetectorNames()) {
                if (geo.getDetector(name).type == "caendt5742") dut_names_.push_back(name);
            }
        }

        for (size_t i = 0; i < dut_names_.size(); ++i) dut_index_[dut_names_[i]] = i;

        // Folder hierarchy: one top folder per configured CFD value, each
        // with a run-wide "total", plus one subfolder per DUT/TREF detector
        // (its own "total" across all channels, plus one "channel_<n>" per
        // physical pixel) - matches corry's per-channel histogram layout.
        // Every histogram's pointer is cached into cfd_total_refs_/
        // dut_total_refs_/dut_channel_refs_ right after registration, so
        // run()'s hot loop never has to rebuild these directory strings or
        // look them up by name again (see ChannelPlotRefs in the header).
        size_t n_cfd = cfd_values_.size();
        size_t n_dut = dut_names_.size();
        cfd_total_refs_.resize(n_cfd);
        dut_total_refs_.assign(n_dut, std::vector<ChannelPlotRefs>(n_cfd));
        dut_channel_refs_.assign(n_dut, std::vector<std::vector<ChannelPlotRefs>>(n_cfd));
        dut_channel_titled_.assign(n_dut, std::vector<std::vector<uint8_t>>(n_cfd));

        for (size_t cfd_idx = 0; cfd_idx < n_cfd; ++cfd_idx) {
            std::string cfd_dir = getName() + "/" + cfdDirName(cfd_values_[cfd_idx]);
            registerChannelPlots(cfd_dir + "/total");
            cfd_total_refs_[cfd_idx] = fetchChannelPlots(cfd_dir + "/total");

            for (size_t dut_idx = 0; dut_idx < n_dut; ++dut_idx) {
                const std::string& det_name = dut_names_[dut_idx];
                if (!geo.hasDetector(det_name)) continue;
                auto& det = geo.getDetector(det_name);

                std::string det_dir = cfd_dir + "/" + det_name;
                registerChannelPlots(det_dir + "/total");
                dut_total_refs_[dut_idx][cfd_idx] = fetchChannelPlots(det_dir + "/total");

                size_t n_channels = static_cast<size_t>(det.n_pixels_x) * static_cast<size_t>(det.n_pixels_y);
                dut_channel_refs_[dut_idx][cfd_idx].resize(n_channels);
                dut_channel_titled_[dut_idx][cfd_idx].assign(n_channels, 0);
                for (int row = 0; row < det.n_pixels_y; ++row) {
                    for (int col = 0; col < det.n_pixels_x; ++col) {
                        int channel = col + det.n_pixels_x * row;
                        std::string channel_dir = det_dir + "/channel_" + std::to_string(channel);
                        registerChannelPlots(channel_dir);
                        auto refs = fetchChannelPlots(channel_dir);
                        // Placeholder until the first real waveform tells us
                        // the raw hardware channel (see #titleChannelPlots's
                        // docs) - still an improvement over the bare
                        // "channel_<n>" name alone, and self-corrects on
                        // first fill if this run never actually sees data
                        // routed here.
                        titleChannelPlots(refs, det_name + " channel_" + std::to_string(channel) +
                                                     " col" + std::to_string(col) + " row" + std::to_string(row) +
                                                     " (raw CH pending)");
                        dut_channel_refs_[dut_idx][cfd_idx][static_cast<size_t>(channel)] = refs;
                    }
                }
            }
        }

        if (save_enabled_) {
            save_h5_ = std::make_unique<H5::H5File>(storage_file_, H5F_ACC_TRUNC);
            save_h5_->createGroup("/waveforms");
            WR_LOG(STATUS, "WaveformProcessingCAEN: saving waveform parameters to " + storage_file_);
        }
    }

    WaveformProcessingCAEN::ChannelWriter& WaveformProcessingCAEN::getChannelWriter(const std::string& det_name) {
        auto it = channel_writers_.find(det_name);
        if (it != channel_writers_.end()) return it->second;

        H5::Group grp = save_h5_->createGroup("/waveforms/" + det_name);
        ChannelWriter w;
        w.event_id   = std::make_unique<AppendableDataset>(grp, "event_id", H5::PredType::NATIVE_UINT64);
        w.col        = std::make_unique<AppendableDataset>(grp, "col", H5::PredType::NATIVE_INT32);
        w.row        = std::make_unique<AppendableDataset>(grp, "row", H5::PredType::NATIVE_INT32);
        w.channel    = std::make_unique<AppendableDataset>(grp, "channel", H5::PredType::NATIVE_INT32);
        w.timestamp  = std::make_unique<AppendableDataset>(grp, "timestamp", H5::PredType::NATIVE_DOUBLE);
        w.baseline   = std::make_unique<AppendableDataset>(grp, "baseline", H5::PredType::NATIVE_DOUBLE);
        w.noise      = std::make_unique<AppendableDataset>(grp, "noise", H5::PredType::NATIVE_DOUBLE);
        w.amplitude  = std::make_unique<AppendableDataset>(grp, "amplitude", H5::PredType::NATIVE_DOUBLE);
        w.snr        = std::make_unique<AppendableDataset>(grp, "snr", H5::PredType::NATIVE_DOUBLE);
        w.rise_time  = std::make_unique<AppendableDataset>(grp, "rise_time", H5::PredType::NATIVE_DOUBLE);
        w.start_time_cfd = std::make_unique<AppendableDataset>(grp, "start_time_cfd", H5::PredType::NATIVE_DOUBLE,
                                                                 /*rank=*/2, /*width=*/cfd_values_.size());
        w.integral   = std::make_unique<AppendableDataset>(grp, "integral", H5::PredType::NATIVE_DOUBLE);
        w.t0         = std::make_unique<AppendableDataset>(grp, "t0", H5::PredType::NATIVE_DOUBLE);
        w.dt         = std::make_unique<AppendableDataset>(grp, "dt", H5::PredType::NATIVE_DOUBLE);
        w.time_over_threshold = std::make_unique<AppendableDataset>(grp, "time_over_threshold", H5::PredType::NATIVE_DOUBLE);
        w.time_over_noise     = std::make_unique<AppendableDataset>(grp, "time_over_noise", H5::PredType::NATIVE_DOUBLE);

        // Records which CFD fraction each start_time_cfd column is, so a
        // reader doesn't have to assume this module's config didn't change
        // between the run that wrote the file and the one reading it back.
        hsize_t n_cfd = cfd_values_.size();
        H5::DataSpace attr_space(1, &n_cfd);
        H5::Attribute attr = grp.createAttribute("cfd_fractions", H5::PredType::NATIVE_DOUBLE, attr_space);
        attr.write(H5::PredType::NATIVE_DOUBLE, cfd_values_.data());

        auto res = channel_writers_.emplace(det_name, std::move(w));
        return res.first->second;
    }

    void WaveformProcessingCAEN::run(DataBatch& batch) {
        auto& geo = GeometryManager::getInstance();
        for (auto& entry : batch.waveforms) {
            const std::string& det_name = entry.first;
            auto& waveforms = entry.second;
            if (waveforms.empty()) continue;

            // Resolved once per detector (not per waveform): looks up this
            // detector's slot into dut_total_refs_/dut_channel_refs_, or
            // SIZE_MAX if det_name isn't a configured DUT.
            size_t dut_idx = SIZE_MAX;
            {
                auto dit = dut_index_.find(det_name);
                if (dit != dut_index_.end()) dut_idx = dit->second;
            }
            int n_pixels_x = geo.hasDetector(det_name) ? geo.getDetector(det_name).n_pixels_x : 0;

            size_t n_wfs = waveforms.size();
            size_t chunk = (n_wfs + thread_pool->getThreadCount() - 1) / thread_pool->getThreadCount();
            size_t n_chunks = (n_wfs + chunk - 1) / chunk;
            std::vector<std::future<void>> futures;

            // Per-thread-local accumulation, merged into one append per
            // detector per batch after all chunks finish - AppendableDataset
            // isn't safe for concurrent writes from multiple threads, same
            // pattern already used for Correlations' parallel histogram
            // fills.
            std::vector<std::vector<uint64_t>> local_event_id(n_chunks);
            std::vector<std::vector<int32_t>> local_col(n_chunks), local_row(n_chunks), local_channel(n_chunks);
            std::vector<std::vector<double>> local_ts(n_chunks), local_baseline(n_chunks), local_noise(n_chunks),
                local_amp(n_chunks), local_snr(n_chunks), local_rise(n_chunks),
                local_integral(n_chunks), local_t0(n_chunks), local_dt(n_chunks),
                local_tot(n_chunks), local_ton(n_chunks);
            // Flattened row-major, width cfd_values_.size() per waveform -
            // one start_time per configured CFD fraction (local_start above
            // keeps just the canonical cfd_values_[0] value, used by the
            // raw_signal_start_time histogram and is_valid_hit gating).
            std::vector<std::vector<double>> local_start_cfd(n_chunks);
            std::vector<double> local_total_time(n_chunks, 0.0), local_fill_time(n_chunks, 0.0);

            size_t chunk_idx = 0;
            for (size_t i = 0; i < n_wfs; i += chunk, ++chunk_idx) {
                size_t end = std::min(i + chunk, n_wfs);
                futures.push_back(thread_pool->submit([&, i, end, chunk_idx]() {
                    auto& pm = PlotManager::getInstance();

                    // Reused across every waveform in this chunk instead of
                    // allocated fresh per waveform - cfd_values_.size() never
                    // changes mid-run, so sizing these once per chunk task
                    // and just overwriting their contents each iteration
                    // avoids a heap allocation per waveform. uint8_t, not
                    // vector<bool>: the bit-packed specialization's per-access
                    // overhead defeats the point of reusing the buffer.
                    // Last two entries are always the dedicated 50%-of-
                    // amplitude and noise-threshold rising-edge crossings
                    // (independent of the configured cfd_values_ list) - the
                    // rising-edge half of time_over_threshold/time_over_noise
                    // (see below); ToT needs a 50% crossing regardless of
                    // whether 0.5 happens to be a configured `cfd` value.
                    size_t tot_idx = 2 + cfd_values_.size();
                    size_t ton_idx = tot_idx + 1;
                    size_t n_finds = ton_idx + 1;
                    std::vector<double> find_thresholds(n_finds), find_results(n_finds);
                    std::vector<uint8_t> find_done(n_finds);
                    // Reused the same way: (voltage, time) samples collected
                    // for the current rising edge's threshold search - see
                    // interpolateTimeAtVoltage() below for what these feed.
                    std::vector<std::pair<double, double>> rising_edge_window, sorted_window;
                    rising_edge_window.reserve(64);
                    sorted_window.reserve(64);
                    // Falling-edge mirror of the above, for the other half of
                    // time_over_threshold/time_over_noise (see the forward
                    // scan after the existing backward one below). Always
                    // exactly 2 thresholds: {tot, ton} - never tied to
                    // cfd_values_, since a falling-edge start_time_cfd has no
                    // meaning (CFD timing is a rising-edge-only concept).
                    std::vector<double> falling_find_thresholds(2), falling_find_results(2);
                    std::vector<uint8_t> falling_find_done(2);
                    std::vector<std::pair<double, double>> falling_edge_window, falling_sorted_window;
                    falling_edge_window.reserve(64);
                    falling_sorted_window.reserve(64);

                    for (size_t w = i; w < end; ++w) {
                        ScopeTimer wf_timer(local_total_time[chunk_idx]);
                        auto& wf = waveforms[w];
                        const auto& data = wf->samples();
                        if (data.size() < 2) continue;

                        Eigen::Map<const Eigen::VectorXd> raw(data.data(), data.size());
                        // ADC code -> Volts (see kCaenAdcMid comment above),
                        // then polarity flip - same net effect as the old
                        // `raw * polarity_` for a non-CAEN waveform source,
                        // since dc_offset_volts defaults to 0.0 and this is
                        // still just an affine transform of `raw` either way.
                        Eigen::VectorXd s = ((raw.array() - kCaenAdcMid) * (kCaenVpp / kCaenAdcLsbDen) + wf->dc_offset_volts).matrix()
                                             * static_cast<double>(polarity_);

                        // Baseline/noise/amplitude, ported from corry's
                        // WfSignal (WfPositiveSignal::_getSignalStartIndex /
                        // _getBaseline / _getNoise / _getAmplitude), rather
                        // than a fixed leading-window estimate. The true
                        // signal onset varies event-to-event and is
                        // typically far from the start of the trace - a
                        // fixed small window both estimates noise from far
                        // fewer samples than corry does (more statistical
                        // scatter in the noise/SNR estimate, which alone can
                        // push otherwise-good pulses outside a tight
                        // two-sided SNR cut) and risks overlapping the actual
                        // rising edge for events where the pulse starts
                        // early, corrupting baseline/amplitude for exactly
                        // those events.
                        Eigen::Index max_idx;
                        s.maxCoeff(&max_idx);
                        if (max_idx < 1) continue; // no samples before the peak at all

                        std::vector<double> before_peak(s.data(), s.data() + max_idx);
                        double median_before_peak = medianOf(before_peak);
                        std::vector<double> deviations(before_peak.size());
                        for (size_t k = 0; k < before_peak.size(); ++k) {
                            deviations[k] = std::abs(before_peak[k] - median_before_peak);
                        }
                        // https://en.wikipedia.org/wiki/Median_absolute_deviation#Relation_to_standard_deviation
                        double mad_before_peak = medianOf(deviations) * 1.4826;

                        // Signal start = the LAST pre-peak sample still
                        // consistent with the pre-peak noise distribution
                        // (median + 1 MAD-sigma) - i.e. the sample right
                        // before the rising edge actually begins.
                        Eigen::Index signal_start_index = -1;
                        double pre_signal_threshold = median_before_peak + mad_before_peak;
                        for (Eigen::Index k = 0; k < max_idx; ++k) {
                            if (s[k] < pre_signal_threshold) signal_start_index = k;
                        }
                        if (signal_start_index < 0) continue; // no valid pre-signal region found

                        // corry's exact off-by-one: baseline/noise are taken
                        // over samples strictly before (signal_start_index -
                        // 1), leaving a 1-sample buffer ahead of the
                        // detected transition.
                        Eigen::Index base_end = signal_start_index - 1;
                        if (base_end <= 0) continue;

                        // Sequential accumulation (not Eigen's .mean()/
                        // .square().mean(), which may use SIMD/pairwise
                        // reduction with different rounding) - matches
                        // corry's std::accumulate order exactly, in case
                        // last-bit differences here get amplified by the
                        // threshold-crossing searches below.
                        double baseline = 0.0;
                        for (Eigen::Index k = 0; k < base_end; ++k) baseline += s[k];
                        baseline /= static_cast<double>(base_end);
                        double noise_var = 0.0;
                        for (Eigen::Index k = 0; k < base_end; ++k) {
                            double d = s[k] - baseline;
                            noise_var += d * d;
                        }
                        noise_var /= static_cast<double>(base_end);
                        double noise = std::sqrt(noise_var);
                        if (noise <= 0.0) continue;

                        // max(sample - baseline) over the WHOLE trace, clamped
                        // at >= 0 (matches corry's loop, which starts
                        // amplitude at 0.0 and only raises it).
                        double amplitude = std::max(0.0, (s.array() - baseline).maxCoeff());
                        double snr = amplitude / noise;

                        wf->baseline = baseline;
                        wf->noise = noise;
                        wf->amplitude = amplitude;
                        wf->snr = snr;

                        // Every threshold crossing needed below (rise_time's
                        // two, plus one per configured CFD value) matches
                        // corry's own findTimeAtRisingEdge(): not a single
                        // 2-point interpolation at the first sample past the
                        // threshold, but a linear interpolation over EVERY
                        // sample from where the pulse leaves its 99%-of-
                        // amplitude peak down to the threshold itself,
                        // sorted by voltage (corry builds a ROOT
                        // TGraph(voltage, time) over that same window and
                        // calls Eval() - see interpolateTimeAtVoltage()'s
                        // docs above for why that reduces to this without
                        // needing ROOT). A naive single-crossing interpolation
                        // is not equivalent: rise_time and snr are both
                        // differences/ratios of two crossings, so a small
                        // per-crossing bias gets amplified in those
                        // quantities even though it barely affects
                        // signalStartTime alone.
                        //
                        // One shared backward walk from the peak builds this
                        // window ONCE and reuses it for every threshold that
                        // hasn't fired yet, rather than corry's own
                        // O(n_thresholds) independent rescans of what's
                        // largely the same overlapping region - every
                        // threshold's own window is exactly what corry's
                        // separate rescan would have collected too (all
                        // windows share the same 99% upper boundary and
                        // start of collection; a lower threshold's window is
                        // simply whatever a higher one's already was, plus
                        // more samples), just computed once.
                        int find_start = std::min(static_cast<int>(max_idx), static_cast<int>(s.size()) - 2);
                        find_thresholds[0] = baseline + 0.9 * amplitude;
                        find_thresholds[1] = baseline + 0.1 * amplitude;
                        for (size_t k = 0; k < cfd_values_.size(); ++k) {
                            find_thresholds[2 + k] = baseline + cfd_values_[k] * amplitude;
                        }
                        // Time-over-threshold's 50%-of-amplitude level, and
                        // time-over-noise's noise level (baseline + N sigma,
                        // N = noise_factor_ - the same knob already used
                        // elsewhere in this function, not the /2-halved
                        // variant used for the charge-integral window below,
                        // which is a different, integral-specific threshold).
                        double tot_thresh_abs = baseline + 0.5 * amplitude;
                        double ton_thresh_abs = baseline + noise_factor_ * noise;
                        find_thresholds[tot_idx] = tot_thresh_abs;
                        find_thresholds[ton_idx] = ton_thresh_abs;
                        std::fill(find_results.begin(), find_results.end(), 0.0);
                        std::fill(find_done.begin(), find_done.end(), 0);
                        size_t n_remaining = n_finds;

                        double peak_thresh_99 = baseline + 0.99 * amplitude;
                        rising_edge_window.clear();
                        // corry's window boundary is off-by-one from what it
                        // looks like at first read: k_stop_rise is set to
                        // k+1 on every k that's STILL above 99%, so once the
                        // loop finally drops below 99% and stops updating
                        // it, k_stop_rise-1 - the largest index actually
                        // included in [start_rise, stop_rise) - is the LAST
                        // sample that was still above 99%, not the first one
                        // at/below it. Tracked here as the single carried-
                        // over point so every threshold's window gets a
                        // genuine bracket even when the very first at/below-
                        // 99% sample already undercuts it (steep pulses that
                        // cross several % thresholds within one sample -
                        // without this point, that collapses to a 1-point
                        // window and can't distinguish between different
                        // thresholds' crossing times at all).
                        bool have_last_above99 = false;
                        double last_above99_v = 0.0, last_above99_t = 0.0;

                        for (int j = find_start; j > 0 && n_remaining > 0; --j) {
                            double v = s[j];
                            if (v > peak_thresh_99) {
                                last_above99_v = v;
                                last_above99_t = j * sampling_interval_;
                                have_last_above99 = true;
                                continue;
                            }
                            if (rising_edge_window.empty() && have_last_above99) {
                                rising_edge_window.emplace_back(last_above99_v, last_above99_t);
                            }
                            rising_edge_window.emplace_back(v, j * sampling_interval_);

                            for (size_t k = 0; k < n_finds; ++k) {
                                if (find_done[k]) continue;
                                if (v < find_thresholds[k]) {
                                    sorted_window = rising_edge_window;
                                    std::sort(sorted_window.begin(), sorted_window.end());
                                    find_results[k] = interpolateTimeAtVoltage(sorted_window, find_thresholds[k]);
                                    find_done[k] = true;
                                    --n_remaining;
                                }
                            }
                        }

                        wf->rise_time = find_results[0] - find_results[1];

                        // Falling-edge mirror of the rising-edge search just
                        // above (Warlock-only addition, no corry equivalent -
                        // corry never computes a falling-edge crossing at
                        // all): walks FORWARD from the peak instead of
                        // backward, using the exact same "carry the last
                        // above-99% point into the window" technique so the
                        // very first post-peak samples (still near 100%,
                        // descending slowly) don't collapse the window to a
                        // single point. Combined with the rising-edge crossings
                        // already found above (find_results[tot_idx]/[ton_idx]),
                        // this gives both edges needed for
                        // time_over_threshold/time_over_noise below.
                        std::fill(falling_find_results.begin(), falling_find_results.end(), 0.0);
                        std::fill(falling_find_done.begin(), falling_find_done.end(), 0);
                        falling_find_thresholds[0] = tot_thresh_abs;
                        falling_find_thresholds[1] = ton_thresh_abs;
                        size_t falling_n_remaining = 2;

                        falling_edge_window.clear();
                        bool have_last_above99_falling = false;
                        double last_above99_v_falling = 0.0, last_above99_t_falling = 0.0;

                        for (int j = static_cast<int>(max_idx); j < s.size() && falling_n_remaining > 0; ++j) {
                            double v = s[j];
                            if (v > peak_thresh_99) {
                                last_above99_v_falling = v;
                                last_above99_t_falling = j * sampling_interval_;
                                have_last_above99_falling = true;
                                continue;
                            }
                            if (falling_edge_window.empty() && have_last_above99_falling) {
                                falling_edge_window.emplace_back(last_above99_v_falling, last_above99_t_falling);
                            }
                            falling_edge_window.emplace_back(v, j * sampling_interval_);

                            for (size_t k = 0; k < 2; ++k) {
                                if (falling_find_done[k]) continue;
                                if (v < falling_find_thresholds[k]) {
                                    falling_sorted_window = falling_edge_window;
                                    std::sort(falling_sorted_window.begin(), falling_sorted_window.end());
                                    falling_find_results[k] = interpolateTimeAtVoltage(falling_sorted_window, falling_find_thresholds[k]);
                                    falling_find_done[k] = true;
                                    --falling_n_remaining;
                                }
                            }
                        }

                        // A pulse that never falls back below threshold
                        // before the trace ends (falling_find_done[k] still
                        // false) leaves falling_find_results[k] at 0.0 from
                        // the fill() above - std::max(..., 0.0) below turns
                        // that into a 0-length ToT/ToN rather than a bogus
                        // negative one, same defensive clamp as corry's own
                        // "no valid signal -> 0, not NaN" convention used for
                        // #integral just below.
                        wf->time_over_threshold = std::max(0.0, falling_find_results[0] - find_results[tot_idx]);
                        wf->time_over_noise = std::max(0.0, falling_find_results[1] - find_results[ton_idx]);

                        // Charge, ported from corry's
                        // WfPositiveSignal::_getSignalIntegralDt()
                        // (WfSignal.cpp) - NOT _getSignalIntegral(): corry's
                        // WaveformProcessingCAEN::run() picks between the two
                        // based on whether `integral_time` is negative
                        // (config default: 2.8, i.e. not negative), so
                        // getSignalIntegralDt() is the one that runs by
                        // default.
                        //
                        // Unlike _getSignalIntegral() (bounded by BOTH a
                        // rising and a falling edge), this one only finds
                        // the rising edge - scanning backward from the peak
                        // for the low threshold (`noise_factor*noise/2`
                        // above baseline; corry expresses it as a fraction
                        // of amplitude but amplitude cancels once converted
                        // to an absolute level) and the 99%-of-amplitude
                        // high threshold - then integrates a FIXED-WIDTH
                        // window of `integral_time_` ns starting there,
                        // regardless of where the pulse actually falls back
                        // down. The window's RAW (baseline-INCLUDED, not
                        // subtracted) samples are integrated with ROOT's
                        // TGraph::Integral() (shoelace formula, closed by a
                        // segment back from the last point to the first,
                        // then |.|*0.5 - see TGraph.cxx), reduced here to
                        // its uniform-time-spacing form instead of building
                        // an actual TGraph.
                        double integral_thresh_abs = baseline + noise_factor_ * noise / 2.0;
                        double rise_high_abs = baseline + amplitude * 0.99;

                        Eigen::Index k_start_rise = -1, k_stop_rise = -1;
                        for (Eigen::Index k = max_idx; k > 0; --k) {
                            double v = s[k];
                            if (v > rise_high_abs) k_stop_rise = k + 1;
                            if (v < integral_thresh_abs) { k_start_rise = k; break; }
                        }

                        // No valid rising edge -> corry returns NaN for this
                        // waveform's signal_integral; left at 0.0 here
                        // instead so it doesn't poison downstream histogram
                        // Mean/StdDev sums (this only fires for the same
                        // edge cases already excluded above via the
                        // noise<=0.0/signal_start_index checks).
                        wf->integral = 0.0;
                        if (k_start_rise >= 0 && k_stop_rise >= 0 && k_start_rise < k_stop_rise) {
                            Eigen::Index first_idx = k_start_rise;
                            // last_integral_index = (integral_time + first_idx*dt)/dt,
                            // truncated via integer cast exactly like corry's
                            // own size_t(...) cast.
                            Eigen::Index last_integral_index = static_cast<Eigen::Index>(
                                (integral_time_ + static_cast<double>(first_idx) * sampling_interval_) / sampling_interval_);
                            Eigen::Index n = last_integral_index - first_idx; // corry's `indices.size()`
                            if (n >= 2 && first_idx + n - 1 < s.size()) {
                                Eigen::Index last_pt = first_idx + n - 1;
                                double sum = 0.0;
                                for (Eigen::Index k = first_idx; k < last_pt; ++k) {
                                    sum += (s[k] + s[k + 1]) * sampling_interval_;
                                }
                                sum += (s[last_pt] + s[first_idx]) * (-static_cast<double>(n - 1) * sampling_interval_);
                                wf->integral = 0.5 * std::abs(sum);
                            }
                        }

                        // Ported from corry's is_valid_hit_(): a waveform
                        // whose snr/start_time/rise_time/integral fall
                        // outside the configured ranges is dropped entirely
                        // (this is what corry's default thresholds actually
                        // gate on in practice - not a real cut, just
                        // rejecting waveforms where no valid signal was
                        // found at all, e.g. integral stuck at 0.0 above
                        // because no rising edge crossed the noise
                        // threshold). Uses cfd_values_[0]'s start_time as
                        // the reference, matching corry's single `cfd`
                        // config value (getSignalStartTime() and
                        // getSignalStartTime_cfd() use the same threshold
                        // there); this only matters in degenerate cases
                        // since min/max_start_time default to (0, huge)
                        // regardless of which cfd fraction is used.
                        bool valid_hit = std::fabs(wf->snr) > min_snr_ && std::fabs(wf->snr) < max_snr_ &&
                                          find_results[2] > min_start_time_ && find_results[2] < max_start_time_ &&
                                          wf->rise_time > min_rise_time_ && wf->rise_time < max_rise_time_ &&
                                          wf->integral > min_integral_ && wf->integral < max_integral_;
                        if (!valid_hit) {
                            continue;
                        }

                        wf->is_processed = true;

                        int channel = n_pixels_x > 0 ? (wf->column() + n_pixels_x * wf->row()) : 0;

                        // start_time is the only CFD-dependent quantity, so
                        // it gets recomputed once per configured CFD value;
                        // the first entry is treated as the canonical value
                        // stored on the waveform itself (used by the flat
                        // raw_* histograms below and by the save file). The
                        // other 6 quantities are identical across every CFD
                        // value for a given waveform - they're re-filled into
                        // each CFD folder's own copy of these histograms
                        // (matching corry's per-CFD directory layout) via
                        // cached Histogram1D pointers (see ChannelPlotRefs),
                        // not by rebuilding "<dir>/<name>" strings and doing
                        // a PlotManager lookup per histogram per CFD value.
                        ScopeTimer fill_timer(local_fill_time[chunk_idx]);
                        wf->start_time_by_cfd.reserve(cfd_values_.size());
                        for (size_t cfd_idx = 0; cfd_idx < cfd_values_.size(); ++cfd_idx) {
                            double start_time = find_results[2 + cfd_idx];
                            if (cfd_idx == 0) wf->start_time = start_time;
                            wf->start_time_by_cfd.emplace_back(cfd_values_[cfd_idx], start_time);
                            local_start_cfd[chunk_idx].push_back(start_time);

                            fillChannelPlots(cfd_total_refs_[cfd_idx], wf->snr, wf->amplitude, wf->noise,
                                              wf->baseline, wf->rise_time, start_time, wf->integral,
                                              wf->time_over_threshold, wf->time_over_noise);

                            if (dut_idx != SIZE_MAX) {
                                fillChannelPlots(dut_total_refs_[dut_idx][cfd_idx], wf->snr, wf->amplitude,
                                                  wf->noise, wf->baseline, wf->rise_time, start_time, wf->integral,
                                                  wf->time_over_threshold, wf->time_over_noise);
                                auto& channel_refs = dut_channel_refs_[dut_idx][cfd_idx];
                                if (static_cast<size_t>(channel) < channel_refs.size()) {
                                    // Stamp this folder's plots with the raw
                                    // CAEN hardware channel on first fill -
                                    // only known here, from actual data, not
                                    // at registration time (see
                                    // titleChannelPlots()'s docs on why
                                    // channel_<n> alone is misleading: it's
                                    // NOT the same number as the physical
                                    // CH<n> label).
                                    if (!dut_channel_titled_[dut_idx][cfd_idx][static_cast<size_t>(channel)]) {
                                        titleChannelPlots(channel_refs[static_cast<size_t>(channel)],
                                            det_name + " CH" + std::to_string(wf->channel()) +
                                            " col" + std::to_string(wf->column()) + " row" + std::to_string(wf->row()));
                                        dut_channel_titled_[dut_idx][cfd_idx][static_cast<size_t>(channel)] = 1;
                                    }
                                    fillChannelPlots(channel_refs[static_cast<size_t>(channel)], wf->snr,
                                                      wf->amplitude, wf->noise, wf->baseline, wf->rise_time,
                                                      start_time, wf->integral, wf->time_over_threshold,
                                                      wf->time_over_noise);
                                }
                            }
                        }

                        {
                            ScopeTimer fill_timer(local_fill_time[chunk_idx]);
                            pm.fill1D(getName(), "raw_snr", wf->snr);
                            pm.fill1D(getName(), "raw_amplitude", wf->amplitude);
                            pm.fill1D(getName(), "raw_rise_time", wf->rise_time);
                            pm.fill1D(getName(), "raw_baseline", wf->baseline);
                            pm.fill1D(getName(), "raw_noise", wf->noise);
                            pm.fill1D(getName(), "raw_signal_start_time", wf->start_time);
                            pm.fill1D(getName(), "raw_integral", wf->integral);
                            pm.fill1D(getName(), "raw_time_over_threshold", wf->time_over_threshold);
                            pm.fill1D(getName(), "raw_time_over_noise", wf->time_over_noise);
                        }

                        if (save_enabled_) {
                            local_event_id[chunk_idx].push_back(wf->eventID());
                            local_col[chunk_idx].push_back(wf->column());
                            local_row[chunk_idx].push_back(wf->row());
                            local_channel[chunk_idx].push_back(wf->channel());
                            local_ts[chunk_idx].push_back(wf->timestamp());
                            local_baseline[chunk_idx].push_back(wf->baseline);
                            local_noise[chunk_idx].push_back(wf->noise);
                            local_amp[chunk_idx].push_back(wf->amplitude);
                            local_snr[chunk_idx].push_back(wf->snr);
                            local_rise[chunk_idx].push_back(wf->rise_time);
                            local_integral[chunk_idx].push_back(wf->integral);
                            local_t0[chunk_idx].push_back(wf->t0);
                            local_dt[chunk_idx].push_back(wf->dt);
                            local_tot[chunk_idx].push_back(wf->time_over_threshold);
                            local_ton[chunk_idx].push_back(wf->time_over_noise);
                        }
                    }
                }));
            }
            for (auto& f : futures) f.get();
            total_processed_ += n_wfs;

            double batch_total_time = 0.0, batch_fill_time = 0.0;
            for (size_t c = 0; c < n_chunks; ++c) {
                batch_total_time += local_total_time[c];
                batch_fill_time += local_fill_time[c];
            }
            stat_time_fill_s += batch_fill_time;
            stat_time_math_s += (batch_total_time - batch_fill_time);

            if (save_enabled_) {
                auto& cw = getChannelWriter(det_name);
                for (size_t c = 0; c < n_chunks; ++c) {
                    cw.event_id->append(local_event_id[c]);
                    cw.col->append(local_col[c]);
                    cw.row->append(local_row[c]);
                    cw.channel->append(local_channel[c]);
                    cw.timestamp->append(local_ts[c]);
                    cw.baseline->append(local_baseline[c]);
                    cw.noise->append(local_noise[c]);
                    cw.amplitude->append(local_amp[c]);
                    cw.snr->append(local_snr[c]);
                    cw.rise_time->append(local_rise[c]);
                    cw.start_time_cfd->append2D(local_start_cfd[c], cfd_values_.size());
                    cw.integral->append(local_integral[c]);
                    cw.t0->append(local_t0[c]);
                    cw.dt->append(local_dt[c]);
                    cw.time_over_threshold->append(local_tot[c]);
                    cw.time_over_noise->append(local_ton[c]);
                }
            }
        }
    }

    void WaveformProcessingCAEN::finalize() {
        WR_LOG(STATUS, "Waveform Math Engine: Processed " + std::to_string(total_processed_) + " pulses.");
        std::ostringstream phase_stats;
        phase_stats << "WaveformProcessingCAEN CPU-time (seconds, summed across worker threads): "
                    << "math=" << stat_time_math_s
                    << " plot_fill=" << stat_time_fill_s;
        WR_LOG(STATUS, phase_stats.str());
        if (save_enabled_ && save_h5_) {
            save_h5_->close();
            WR_LOG(STATUS, "WaveformProcessingCAEN: saved parameters to " + storage_file_);
        }
    }
    REGISTER_MODULE(WaveformProcessingCAEN)
}