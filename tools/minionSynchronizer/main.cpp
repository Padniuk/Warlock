/**
 * @file main.cpp
 * @brief CLI entry point for minionSynchronizer - see Synchronizer.hpp for
 * the three-stage detection logic (per-detector CUSUM -> shared-onset max
 * -> pooled banded shift search) this drives.
 *
 * Two-phase workflow, deliberately never fully automatic (this edits real
 * experiment data): (1) default mode detects and REPORTS the group's
 * shared onset and pooled shift trajectory, writing nothing; (2) once
 * you've reviewed that, -y applies it directly to every detector in -d's
 * group.
 */
#include "Synchronizer.hpp"
#include <H5Cpp.h>
#include <getopt.h>
#include <iostream>
#include <sstream>

namespace {
    void printUsage(const char* prog) {
        std::cerr << "Usage: " << prog << " <tracks.h5> <waveforms.h5> <geometry.toml> [options]\n"
                  << "\n"
                  << "  Detect mode (default): for -d's group (every detector sharing one\n"
                  << "  physical digitizer), runs CUSUM independently per detector, combines\n"
                  << "  to one shared onset (the MAX of the individually-confirmed break\n"
                  << "  points - a too-early individual trigger is far more likely than a\n"
                  << "  too-late one), then a banded shift search POOLED across the whole\n"
                  << "  group finds the shift trajectory from that onset onward. Writes\n"
                  << "  nothing.\n"
                  << "\n"
                  << "  If -d is omitted, every detector present in both the geometry and\n"
                  << "  waveforms.h5 is checked, but each as its OWN independent group of 1\n"
                  << "  (no cross-detector combination) - pass a comma-separated -d to\n"
                  << "  actually group detectors that share a digitizer.\n"
                  << "\n"
                  << "  Apply mode (-y): shifts waveforms.h5's real event_id column for\n"
                  << "  EVERY detector in -d's group using THIS SAME run's own detected\n"
                  << "  trajectory, and writes a new file - never the input, unless -p is\n"
                  << "  also given.\n"
                  << "\n"
                  << "Required:\n"
                  << "  -d <det1,det2,...>  The detector group to analyze jointly (required for\n"
                  << "                      apply mode; comma-separated for a real group, a\n"
                  << "                      single name for one detector; omit entirely to check\n"
                  << "                      every detector independently)\n"
                  << "\n"
                  << "WaveformSelector-equivalent cuts (same names/semantics/defaults as the\n"
                  << "real module - set these to whatever the real analysis config uses, same\n"
                  << "for every detector in the group):\n"
                  << "  --min-snr <v>        (default 0)      --max-snr <v>        (default off)\n"
                  << "  --min-amplitude <v>  (default 0)      --max-amplitude <v>  (default off)\n"
                  << "  --min-noise <v>      (default 0)      --max-noise <v>      (default off)\n"
                  << "  --min-rise-time <v>  (default 0)      --max-rise-time <v>  (default off)\n"
                  << "  --min-start-time <v> (default 0)      --max-start-time <v> (default off)\n"
                  << "  --min-charge <v>     (default 0)      --max-charge <v>     (default off)\n"
                  << "  --cfd <v>            (default 0.5, must match a saved cfd fraction exactly)\n"
                  << "  --sensoredge <v>     Acceptance tolerance (default 0, DetectorGeo::\n"
                  << "                       hasIntercept()'s own pixel_tolerance)\n"
                  << "\n"
                  << "CUSUM changepoint detection (run per detector, before combination):\n"
                  << "  --cusum-baseline <n> Leading in-acceptance tracks used to estimate the\n"
                  << "                       healthy pass rate (default 5000)\n"
                  << "  --cusum-slack <v>    Ignore deviations smaller than this (default 0.01)\n"
                  << "  --cusum-threshold <v> Cumulative-sum level confirming a persistent drop\n"
                  << "                       (default 15.0)\n"
                  << "\n"
                  << "Pooled banded sequential shift search (from the shared onset onward):\n"
                  << "  --shift-window <n>   Pooled in-acceptance tracks per window, summed\n"
                  << "                       across the group (default 5000)\n"
                  << "  --shift-band <n>     Search +/- this many shifts around the PREVIOUS\n"
                  << "                       window's answer (default 30)\n"
                  << "\n"
                  << "Apply mode:\n"
                  << "  -y                   Use THIS run's own detected trajectory for -d's group\n"
                  << "  -o <output.h5>       Output path (default: <waveforms>_corrected.h5)\n"
                  << "  -p                   Allow overwriting waveforms.h5 in place (refused by\n"
                  << "                       default - a real, hard-to-reverse edit)\n";
    }

    enum LongOpt {
        OPT_MIN_SNR = 1000, OPT_MAX_SNR, OPT_MIN_AMPLITUDE, OPT_MAX_AMPLITUDE,
        OPT_MIN_NOISE, OPT_MAX_NOISE, OPT_MIN_RISE_TIME, OPT_MAX_RISE_TIME,
        OPT_MIN_START_TIME, OPT_MAX_START_TIME, OPT_MIN_CHARGE, OPT_MAX_CHARGE,
        OPT_CFD, OPT_SENSOREDGE, OPT_CUSUM_BASELINE, OPT_CUSUM_SLACK, OPT_CUSUM_THRESHOLD,
        OPT_SHIFT_WINDOW, OPT_SHIFT_BAND,
    };

    std::vector<std::string> splitComma(const std::string& s) {
        std::vector<std::string> out;
        std::stringstream ss(s);
        std::string part;
        while (std::getline(ss, part, ',')) {
            if (!part.empty()) out.push_back(part);
        }
        return out;
    }

    void printGroupResult(const std::string& prog, const framework::SynchronizerConfig& cfg,
                           const framework::GroupResult& g) {
        std::cout << "\n=== group: ";
        for (size_t i = 0; i < cfg.detectors.size(); ++i) std::cout << (i ? "," : "") << cfg.detectors[i];
        std::cout << " ===\n";

        std::cout << "  Per-detector CUSUM (independent):\n";
        for (auto const& pdb : g.per_detector) {
            if (!pdb.shift_tested) {
                std::cout << "    " << pdb.detector << ": no usable tracks/waveforms - skipped\n";
            } else if (!pdb.found_break) {
                std::cout << "    " << pdb.detector << ": baseline=" << pdb.baseline_pass_rate
                          << ", no persistent drop found (not evidence of clean - could be insufficient statistics)\n";
            } else {
                std::cout << "    " << pdb.detector << ": baseline=" << pdb.baseline_pass_rate
                          << ", break at event " << pdb.break_event << "\n";
            }
        }

        if (!g.found_break) {
            std::cout << "  No member of this group found a break - nothing to fix here.\n";
            return;
        }
        std::cout << "  Shared onset (max over members that found a break): event " << g.shared_break_event << "\n"
                   << "  Pooled shift trajectory (" << g.segments.size() << " window(s)):\n";
        for (auto const& seg : g.segments) {
            std::cout << "    [" << seg.event_lo << "-" << seg.event_hi << "] shift=" << seg.shift
                       << "  recovered=" << seg.recovered_rate << "  unshifted=" << seg.unshifted_rate << "\n";
        }
        std::cout << "  To apply: " << prog << " <tracks.h5> <waveforms.h5> <geometry.toml> -d ";
        for (size_t i = 0; i < cfg.detectors.size(); ++i) std::cout << (i ? "," : "") << cfg.detectors[i];
        std::cout << " -y -o <output.h5> [same cut/cfd flags as above]\n";
    }
}

int main(int argc, char** argv) {
    framework::SynchronizerConfig cfg;
    bool use_detected = false;
    std::string output_file;
    bool allow_in_place = false;
    std::string detector_arg;

    static struct option long_opts[] = {
        {"min-snr", required_argument, nullptr, OPT_MIN_SNR},
        {"max-snr", required_argument, nullptr, OPT_MAX_SNR},
        {"min-amplitude", required_argument, nullptr, OPT_MIN_AMPLITUDE},
        {"max-amplitude", required_argument, nullptr, OPT_MAX_AMPLITUDE},
        {"min-noise", required_argument, nullptr, OPT_MIN_NOISE},
        {"max-noise", required_argument, nullptr, OPT_MAX_NOISE},
        {"min-rise-time", required_argument, nullptr, OPT_MIN_RISE_TIME},
        {"max-rise-time", required_argument, nullptr, OPT_MAX_RISE_TIME},
        {"min-start-time", required_argument, nullptr, OPT_MIN_START_TIME},
        {"max-start-time", required_argument, nullptr, OPT_MAX_START_TIME},
        {"min-charge", required_argument, nullptr, OPT_MIN_CHARGE},
        {"max-charge", required_argument, nullptr, OPT_MAX_CHARGE},
        {"cfd", required_argument, nullptr, OPT_CFD},
        {"sensoredge", required_argument, nullptr, OPT_SENSOREDGE},
        {"cusum-baseline", required_argument, nullptr, OPT_CUSUM_BASELINE},
        {"cusum-slack", required_argument, nullptr, OPT_CUSUM_SLACK},
        {"cusum-threshold", required_argument, nullptr, OPT_CUSUM_THRESHOLD},
        {"shift-window", required_argument, nullptr, OPT_SHIFT_WINDOW},
        {"shift-band", required_argument, nullptr, OPT_SHIFT_BAND},
        {nullptr, 0, nullptr, 0},
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "d:yo:p", long_opts, nullptr)) != -1) {
        switch (opt) {
            case 'd': detector_arg = optarg; break;
            case 'y': use_detected = true; break;
            case 'o': output_file = optarg; break;
            case 'p': allow_in_place = true; break;
            case OPT_MIN_SNR: cfg.min_snr = std::stod(optarg); break;
            case OPT_MAX_SNR: cfg.max_snr = std::stod(optarg); break;
            case OPT_MIN_AMPLITUDE: cfg.min_amplitude = std::stod(optarg); break;
            case OPT_MAX_AMPLITUDE: cfg.max_amplitude = std::stod(optarg); break;
            case OPT_MIN_NOISE: cfg.min_noise = std::stod(optarg); break;
            case OPT_MAX_NOISE: cfg.max_noise = std::stod(optarg); break;
            case OPT_MIN_RISE_TIME: cfg.min_rise_time = std::stod(optarg); break;
            case OPT_MAX_RISE_TIME: cfg.max_rise_time = std::stod(optarg); break;
            case OPT_MIN_START_TIME: cfg.min_start_time = std::stod(optarg); break;
            case OPT_MAX_START_TIME: cfg.max_start_time = std::stod(optarg); break;
            case OPT_MIN_CHARGE: cfg.min_charge = std::stod(optarg); break;
            case OPT_MAX_CHARGE: cfg.max_charge = std::stod(optarg); break;
            case OPT_CFD: cfg.cfd = std::stod(optarg); break;
            case OPT_SENSOREDGE: cfg.spatial_cut_sensoredge = std::stod(optarg); break;
            case OPT_CUSUM_BASELINE: cfg.cusum_baseline_tracks = static_cast<size_t>(std::stoull(optarg)); break;
            case OPT_CUSUM_SLACK: cfg.cusum_slack = std::stod(optarg); break;
            case OPT_CUSUM_THRESHOLD: cfg.cusum_threshold = std::stod(optarg); break;
            case OPT_SHIFT_WINDOW: cfg.shift_window_tracks = static_cast<size_t>(std::stoull(optarg)); break;
            case OPT_SHIFT_BAND: cfg.shift_band = std::stoll(optarg); break;
            default: printUsage(argv[0]); return 1;
        }
    }

    if (argc - optind < 3) {
        std::cerr << "[minion] Error: 3 positional arguments required.\n\n";
        printUsage(argv[0]);
        return 1;
    }
    cfg.tracks_file = argv[optind];
    cfg.waveforms_file = argv[optind + 1];
    cfg.geo_file = argv[optind + 2];

    if (use_detected && detector_arg.empty()) {
        std::cerr << "[minion] Error: apply mode needs -d <detector[,detector,...]> - shifting is per-group.\n";
        return 1;
    }

    try {
        if (!detector_arg.empty()) {
            cfg.detectors = splitComma(detector_arg);
            framework::Synchronizer sync(cfg);
            auto g = sync.detect();

            if (use_detected) {
                if (!g.found_break || g.segments.empty()) {
                    std::cerr << "[minion] Error: -y needs a detected, shift-tested break for this group - "
                                 "none found. Run detect mode first to see why.\n";
                    return 1;
                }
                std::string out = output_file.empty() ? cfg.waveforms_file + "_corrected.h5" : output_file;
                std::cout << "[minion] Applying " << g.segments.size() << "-segment pooled shift trajectory to "
                          << cfg.detectors.size() << " detector(s), starting at shared onset event "
                          << g.shared_break_event << ":\n";
                for (auto const& seg : g.segments) {
                    std::cout << "    [" << seg.event_lo << "-" << seg.event_hi << "] shift=" << seg.shift
                              << " (association " << seg.unshifted_rate << " -> " << seg.recovered_rate << ")\n";
                }
                sync.applyShift(cfg.detectors, g.segments, out, allow_in_place);
                return 0;
            }

            printGroupResult(argv[0], cfg, g);
            return 0;
        }

        // No -d: discover every detector and check each independently
        // (a group of 1 - no cross-detector combination).
        framework::Synchronizer discovery(cfg);
        discovery.loadGeometryForDiscovery();
        auto names = discovery.availableDetectors();
        if (names.empty()) {
            std::cerr << "[minion] No detector present in both " << cfg.geo_file << " and " << cfg.waveforms_file << ".\n";
            return 1;
        }
        for (auto const& name : names) {
            cfg.detectors = {name};
            framework::Synchronizer sync(cfg);
            auto g = sync.detect();
            printGroupResult(argv[0], cfg, g);
        }
        return 0;
    } catch (const H5::Exception& e) {
        // H5::Exception doesn't derive from std::exception in the HDF5
        // C++ API, so it isn't caught by the handler below at all - an
        // uncaught one crashes via std::terminate() instead of printing a
        // clean message. Caught explicitly here as defense in depth.
        std::cerr << "[minion] Fatal HDF5 error: " << e.getDetailMsg() << std::endl;
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "[minion] Fatal error: " << e.what() << std::endl;
        return 1;
    }
}
