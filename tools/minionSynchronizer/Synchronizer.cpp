/** @file Synchronizer.cpp */
#include "Synchronizer.hpp"
#include "core/detector/GeometryManager.hpp"
#include <H5Cpp.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <unordered_set>

namespace framework {

    Synchronizer::Synchronizer(SynchronizerConfig cfg) : cfg_(std::move(cfg)) {}

    std::vector<Synchronizer::TrackRow> Synchronizer::loadTracks() const {
        H5::Exception::dontPrint();
        H5::H5File file(cfg_.tracks_file, H5F_ACC_RDONLY);
        H5::Group grp = file.openGroup("tracks");

        auto readCol = [&](const std::string& name, auto& out) {
            H5::DataSet ds = grp.openDataSet(name);
            hsize_t dim[1];
            ds.getSpace().getSimpleExtentDims(dim);
            out.resize(dim[0]);
            using T = typename std::decay_t<decltype(out)>::value_type;
            ds.read(out.data(), std::is_same_v<T, uint64_t> ? H5::PredType::NATIVE_UINT64 : H5::PredType::NATIVE_DOUBLE);
        };
        std::vector<uint64_t> event_id;
        std::vector<double> mx, cx, my, cy;
        readCol("event_id", event_id);
        readCol("mx", mx);
        readCol("cx", cx);
        readCol("my", my);
        readCol("cy", cy);

        std::vector<TrackRow> rows(event_id.size());
        for (size_t i = 0; i < event_id.size(); ++i) rows[i] = {event_id[i], mx[i], cx[i], my[i], cy[i]};
        std::sort(rows.begin(), rows.end(), [](TrackRow const& a, TrackRow const& b) { return a.event_id < b.event_id; });
        return rows;
    }

    std::vector<Synchronizer::WfRow> Synchronizer::loadWaveforms(const std::string& detector) const {
        H5::Exception::dontPrint();
        H5::H5File file(cfg_.waveforms_file, H5F_ACC_RDONLY);
        std::string path = "waveforms/" + detector;
        if (!file.nameExists(path)) return {};
        H5::Group grp = file.openGroup(path);

        H5::DataSet ds_ev = grp.openDataSet("event_id");
        H5::DataSet ds_snr = grp.openDataSet("snr");
        H5::DataSet ds_amp = grp.openDataSet("amplitude");
        H5::DataSet ds_noise = grp.openDataSet("noise");
        H5::DataSet ds_rise = grp.openDataSet("rise_time");
        H5::DataSet ds_charge = grp.openDataSet("integral");
        H5::DataSet ds_start = grp.openDataSet("start_time_cfd");

        hsize_t dim[1];
        ds_ev.getSpace().getSimpleExtentDims(dim);
        size_t n = dim[0];

        std::vector<uint64_t> event_id(n);
        std::vector<double> snr(n), amplitude(n), noise(n), rise_time(n), charge(n);
        ds_ev.read(event_id.data(), H5::PredType::NATIVE_UINT64);
        ds_snr.read(snr.data(), H5::PredType::NATIVE_DOUBLE);
        ds_amp.read(amplitude.data(), H5::PredType::NATIVE_DOUBLE);
        ds_noise.read(noise.data(), H5::PredType::NATIVE_DOUBLE);
        ds_rise.read(rise_time.data(), H5::PredType::NATIVE_DOUBLE);
        ds_charge.read(charge.data(), H5::PredType::NATIVE_DOUBLE);

        hsize_t dim2[2];
        ds_start.getSpace().getSimpleExtentDims(dim2);
        size_t n_cfd = dim2[1];
        std::vector<double> start_all(dim2[0] * dim2[1]);
        ds_start.read(start_all.data(), H5::PredType::NATIVE_DOUBLE);

        int cfd_idx = -1;
        if (grp.attrExists("cfd_fractions")) {
            H5::Attribute attr = grp.openAttribute("cfd_fractions");
            hsize_t adim[1];
            attr.getSpace().getSimpleExtentDims(adim);
            std::vector<double> fractions(adim[0]);
            attr.read(H5::PredType::NATIVE_DOUBLE, fractions.data());
            for (size_t i = 0; i < fractions.size(); ++i) {
                if (std::abs(fractions[i] - cfg_.cfd) < 1e-9) { cfd_idx = static_cast<int>(i); break; }
            }
        }
        if (cfd_idx < 0) {
            throw std::runtime_error("Synchronizer: cfd=" + std::to_string(cfg_.cfd) +
                                      " not found in '" + detector + "'s saved cfd_fractions.");
        }

        std::vector<WfRow> rows(n);
        for (size_t i = 0; i < n; ++i) {
            rows[i] = {event_id[i], snr[i], amplitude[i], noise[i], rise_time[i],
                       start_all[i * n_cfd + static_cast<size_t>(cfd_idx)], charge[i]};
        }
        return rows;
    }

    bool Synchronizer::passesCuts(const WfRow& w) const {
        return w.snr >= cfg_.min_snr && w.snr <= cfg_.max_snr &&
               w.amplitude >= cfg_.min_amplitude && w.amplitude <= cfg_.max_amplitude &&
               w.noise >= cfg_.min_noise && w.noise <= cfg_.max_noise &&
               w.rise_time >= cfg_.min_rise_time && w.rise_time <= cfg_.max_rise_time &&
               w.start_time >= cfg_.min_start_time && w.start_time <= cfg_.max_start_time &&
               w.charge >= cfg_.min_charge && w.charge <= cfg_.max_charge;
    }

    void Synchronizer::loadGeometryForDiscovery() const {
        GeometryManager::getInstance().loadGeometry(cfg_.geo_file);
    }

    std::vector<std::string> Synchronizer::availableDetectors() const {
        H5::Exception::dontPrint();
        H5::H5File file(cfg_.waveforms_file, H5F_ACC_RDONLY);
        if (!file.nameExists("waveforms")) return {};
        H5::Group root = file.openGroup("waveforms");

        std::vector<std::string> names;
        for (hsize_t i = 0; i < root.getNumObjs(); ++i) {
            std::string name = root.getObjnameByIdx(i);
            if (GeometryManager::getInstance().hasDetector(name)) names.push_back(name);
        }
        return names;
    }

    GroupResult Synchronizer::detect() {
        GeometryManager::getInstance().loadGeometry(cfg_.geo_file);
        auto tracks = loadTracks();

        GroupResult group;

        struct DetData {
            std::vector<uint64_t> acc_events;  // sorted, every in-acceptance track's event_id
            std::unordered_set<uint64_t> good_events;
        };
        std::vector<DetData> det_data(cfg_.detectors.size());

        // --- Stage 1: independent per-detector CUSUM ---
        for (size_t di = 0; di < cfg_.detectors.size(); ++di) {
            const std::string& det_name = cfg_.detectors[di];
            PerDetectorBreak pdb;
            pdb.detector = det_name;

            if (!GeometryManager::getInstance().hasDetector(det_name)) {
                group.per_detector.push_back(pdb);
                continue;
            }
            auto& det = GeometryManager::getInstance().getDetector(det_name);
            double z = det.position.z();

            std::vector<uint64_t> acc_events;
            acc_events.reserve(tracks.size());
            for (auto const& t : tracks) {
                Eigen::Vector3d global_pred(t.mx * z + t.cx, t.my * z + t.cy, z);
                Eigen::Vector3d local_pred = det.R_inv * (global_pred - det.position);
                if (det.hasIntercept(local_pred.x(), local_pred.y(), cfg_.spatial_cut_sensoredge))
                    acc_events.push_back(t.event_id);
            }
            if (acc_events.empty()) { group.per_detector.push_back(pdb); continue; }

            auto wf = loadWaveforms(det_name);
            if (wf.empty()) { group.per_detector.push_back(pdb); continue; }

            std::unordered_set<uint64_t> good_events;
            for (auto const& w : wf) {
                if (passesCuts(w)) good_events.insert(w.event_id);
            }
            pdb.shift_tested = true;

            size_t n = acc_events.size();
            std::vector<double> passed(n);
            for (size_t i = 0; i < n; ++i) passed[i] = good_events.count(acc_events[i]) ? 1.0 : 0.0;

            size_t baseline_n = std::min(cfg_.cusum_baseline_tracks, std::max<size_t>(n / 5, 1));
            double target = 0.0;
            for (size_t i = 0; i < baseline_n; ++i) target += passed[i];
            target /= static_cast<double>(baseline_n);
            pdb.baseline_pass_rate = target;

            double C = 0.0;
            int64_t detect_idx = -1;
            for (size_t i = 0; i < n; ++i) {
                double s = target - passed[i] - cfg_.cusum_slack;
                C = std::max(0.0, C + s);
                if (detect_idx < 0 && C > cfg_.cusum_threshold) detect_idx = static_cast<int64_t>(i);
            }
            if (detect_idx >= 0) {
                pdb.found_break = true;
                pdb.break_event = acc_events[static_cast<size_t>(detect_idx)];
            }

            group.per_detector.push_back(pdb);
            det_data[di].acc_events = std::move(acc_events);
            det_data[di].good_events = std::move(good_events);
        }

        // --- Stage 2: shared onset = max over members that found a break ---
        // A detector that never crossed threshold is excluded from the
        // vote entirely - "never fired" means "not enough evidence," not
        // "confirmed clean," so it must never pull the shared estimate
        // earlier by being silently treated as a 0.
        bool any_break = false;
        uint64_t shared = 0;
        for (auto const& pdb : group.per_detector) {
            if (pdb.found_break) {
                any_break = true;
                shared = std::max(shared, pdb.break_event);
            }
        }
        if (!any_break) return group;
        group.found_break = true;
        group.shared_break_event = shared;

        // --- Stage 3: pooled banded shift search across the whole group ---
        struct Entry {
            uint64_t event_id;
            size_t det_idx;
        };
        std::vector<Entry> combined;
        for (size_t di = 0; di < det_data.size(); ++di) {
            for (uint64_t e : det_data[di].acc_events) {
                if (e >= shared) combined.push_back({e, di});
            }
        }
        std::sort(combined.begin(), combined.end(), [](Entry const& a, Entry const& b) { return a.event_id < b.event_id; });

        size_t n = combined.size();
        size_t idx = 0;
        int64_t shift_center = 0;
        while (idx < n) {
            size_t end = std::min(idx + cfg_.shift_window_tracks, n);

            double unshifted = 0.0;
            for (size_t k = idx; k < end; ++k) {
                auto const& e = combined[k];
                if (det_data[e.det_idx].good_events.count(e.event_id)) unshifted += 1.0;
            }
            unshifted /= static_cast<double>(end - idx);

            int64_t best_shift = shift_center;
            double best_rate = -1.0;
            for (int64_t s = shift_center - cfg_.shift_band; s <= shift_center + cfg_.shift_band; ++s) {
                size_t hits = 0;
                for (size_t k = idx; k < end; ++k) {
                    auto const& e = combined[k];
                    int64_t shifted = static_cast<int64_t>(e.event_id) + s;
                    if (shifted >= 0 && det_data[e.det_idx].good_events.count(static_cast<uint64_t>(shifted))) ++hits;
                }
                double rate = static_cast<double>(hits) / static_cast<double>(end - idx);
                if (rate > best_rate) { best_rate = rate; best_shift = s; }
            }

            ShiftSegment seg;
            seg.event_lo = combined[idx].event_id;
            seg.event_hi = combined[end - 1].event_id;
            seg.shift = best_shift;
            seg.recovered_rate = best_rate;
            seg.unshifted_rate = unshifted;
            group.segments.push_back(seg);

            shift_center = best_shift;
            idx = end;
        }

        return group;
    }

    void Synchronizer::applyShift(const std::vector<std::string>& detectors, const std::vector<ShiftSegment>& segments,
                                   const std::string& output_file, bool allow_in_place) {
        if (output_file == cfg_.waveforms_file && !allow_in_place) {
            throw std::runtime_error("Refusing to overwrite " + cfg_.waveforms_file +
                                      " in place without allow_in_place=true - pass a different output "
                                      "path, or explicitly request in-place.");
        }
        if (segments.empty()) {
            throw std::runtime_error("Synchronizer::applyShift: empty shift schedule - nothing to apply.");
        }

        H5::Exception::dontPrint();

        // Always write through a temporary path, then move it into place
        // at the very end - HDF5 can't truncate-and-recreate a file that's
        // still open for reading in the same process, which is exactly
        // what opening `output_file` directly here would do whenever it's
        // the same path as the input (the -p in-place case). Writing to a
        // separate temp file and renaming it over output_file once
        // everything below has succeeded also means a failure partway
        // through never leaves output_file half-written - the original
        // (if output_file == the input) or any prior output stays intact
        // until the rename, which is atomic.
        std::string temp_file = output_file + ".partial";

        {
        H5::H5File in(cfg_.waveforms_file, H5F_ACC_RDONLY);
        H5::H5File out(temp_file, H5F_ACC_TRUNC);

        if (!in.nameExists("waveforms")) throw std::runtime_error("No 'waveforms' group in " + cfg_.waveforms_file);
        H5::Group in_root = in.openGroup("waveforms");
        H5::Group out_root = out.createGroup("waveforms");

        for (hsize_t i = 0; i < in_root.getNumObjs(); ++i) {
            std::string det_name = in_root.getObjnameByIdx(i);
            bool is_target = std::find(detectors.begin(), detectors.end(), det_name) != detectors.end();
            H5::Group in_det = in_root.openGroup(det_name);
            H5::Group out_det = out_root.createGroup(det_name);

            std::vector<uint64_t> event_id;
            {
                H5::DataSet ds = in_det.openDataSet("event_id");
                hsize_t dim[1];
                ds.getSpace().getSimpleExtentDims(dim);
                event_id.resize(dim[0]);
                ds.read(event_id.data(), H5::PredType::NATIVE_UINT64);
            }
            size_t shifted_rows = 0;
            if (is_target) {
                for (auto& e : event_id) {
                    for (auto const& seg : segments) {
                        if (e < seg.event_lo || e > seg.event_hi) continue;
                        // detect()'s search finds shift s such that
                        // track_id + s == waveform_id (see Stage 3 above), i.e.
                        // waveform_id is s ahead of track_id. To make the
                        // waveform's OWN stored event_id line up with the
                        // track directly (so a later shift=0 lookup matches),
                        // this file's event_id must move by -s, not +s -
                        // applying +s here doubles the mismatch instead of
                        // canceling it.
                        int64_t shifted = static_cast<int64_t>(e) - seg.shift;
                        if (shifted < 0) {
                            throw std::runtime_error("Shift would make event_id negative for a row in '" +
                                                      det_name + "' - refusing to write a corrupt file.");
                        }
                        e = static_cast<uint64_t>(shifted);
                        ++shifted_rows;
                        break;  // segments cover disjoint event_id ranges - first match wins
                    }
                }
                std::cout << "[minion] '" << det_name << "': shifted " << shifted_rows << " row(s) across "
                          << segments.size() << " segment(s)." << std::endl;
            }
            {
                hsize_t dim[1] = {event_id.size()};
                H5::DataSpace space(1, dim);
                H5::DataSet ds = out_det.createDataSet("event_id", H5::PredType::NATIVE_UINT64, space);
                ds.write(event_id.data(), H5::PredType::NATIVE_UINT64);
            }

            static const std::vector<std::string> int32_fields = {"col", "row", "channel"};
            static const std::vector<std::string> double_fields = {"timestamp", "baseline", "noise", "amplitude",
                                                                      "snr", "rise_time", "integral", "t0", "dt"};
            for (auto const& name : int32_fields) {
                if (!in_det.nameExists(name)) continue;
                H5::DataSet ds_in = in_det.openDataSet(name);
                hsize_t dim[1];
                ds_in.getSpace().getSimpleExtentDims(dim);
                std::vector<int32_t> buf(dim[0]);
                ds_in.read(buf.data(), H5::PredType::NATIVE_INT32);
                H5::DataSpace space(1, dim);
                H5::DataSet ds_out = out_det.createDataSet(name, H5::PredType::NATIVE_INT32, space);
                ds_out.write(buf.data(), H5::PredType::NATIVE_INT32);
            }
            for (auto const& name : double_fields) {
                if (!in_det.nameExists(name)) continue;
                H5::DataSet ds_in = in_det.openDataSet(name);
                hsize_t dim[1];
                ds_in.getSpace().getSimpleExtentDims(dim);
                std::vector<double> buf(dim[0]);
                ds_in.read(buf.data(), H5::PredType::NATIVE_DOUBLE);
                H5::DataSpace space(1, dim);
                H5::DataSet ds_out = out_det.createDataSet(name, H5::PredType::NATIVE_DOUBLE, space);
                ds_out.write(buf.data(), H5::PredType::NATIVE_DOUBLE);
            }
            if (in_det.nameExists("start_time_cfd")) {
                H5::DataSet ds_in = in_det.openDataSet("start_time_cfd");
                hsize_t dim[2];
                ds_in.getSpace().getSimpleExtentDims(dim);
                std::vector<double> buf(dim[0] * dim[1]);
                ds_in.read(buf.data(), H5::PredType::NATIVE_DOUBLE);
                H5::DataSpace space(2, dim);
                H5::DataSet ds_out = out_det.createDataSet("start_time_cfd", H5::PredType::NATIVE_DOUBLE, space);
                ds_out.write(buf.data(), H5::PredType::NATIVE_DOUBLE);
            }
            if (in_det.attrExists("cfd_fractions")) {
                H5::Attribute attr_in = in_det.openAttribute("cfd_fractions");
                hsize_t adim[1];
                attr_in.getSpace().getSimpleExtentDims(adim);
                std::vector<double> fractions(adim[0]);
                attr_in.read(H5::PredType::NATIVE_DOUBLE, fractions.data());
                H5::DataSpace aspace(1, adim);
                H5::Attribute attr_out = out_det.createAttribute("cfd_fractions", H5::PredType::NATIVE_DOUBLE, aspace);
                attr_out.write(H5::PredType::NATIVE_DOUBLE, fractions.data());
            }
        }
        }  // `in`/`out` close here (RAII) - required before the rename below.

        if (std::rename(temp_file.c_str(), output_file.c_str()) != 0) {
            throw std::runtime_error("Failed to move " + temp_file + " to " + output_file + ": " + std::strerror(errno));
        }

        std::cout << "[minion] Wrote " << output_file << " for " << detectors.size() << " detector(s)." << std::endl;
    }
}
