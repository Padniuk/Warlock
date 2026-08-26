#include "core/detector/GeometryManager.hpp"
#include "core/detector/layouts/OrthogonalLayout.hpp"
#include "core/config/Configuration.hpp"
#include "core/utils/Units.hpp"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <cmath>

namespace framework {

    namespace {
        // Winding-number point-in-polygon test, ported verbatim from corry's
        // PixelDetector::winding_number/isLeft (http://geomalgorithms.com/a03-_inclusion.html).
        int isLeft(std::pair<int, int> p0, std::pair<int, int> p1, std::pair<int, int> p2) {
            return (p1.first - p0.first) * (p2.second - p0.second) - (p2.first - p0.first) * (p1.second - p0.second);
        }
        int windingNumber(std::pair<int, int> probe, const std::vector<std::pair<int, int>>& polygon) {
            if (polygon.size() < 3) return 0;
            int wn = 0;
            for (size_t i = 0; i < polygon.size(); ++i) {
                auto p0 = polygon[i];
                auto p1 = (i + 1 < polygon.size()) ? polygon[i + 1] : polygon[0];
                if (p0.second <= probe.second) {
                    if (p1.second > probe.second && isLeft(p0, p1, probe) > 0) ++wn;
                } else {
                    if (p1.second <= probe.second && isLeft(p0, p1, probe) < 0) --wn;
                }
            }
            return wn;
        }
    }

    bool DetectorGeo::hasIntercept(double local_x, double local_y, double pixel_tolerance) const {
        double column = getColumn(local_x);
        if (column < pixel_tolerance - 0.5 || column > (n_pixels_x - pixel_tolerance - 0.5)) {
            return false;
        }
        // Y can't reuse the same "row units" check as X once a multi-die
        // board has a real inter-die gap: getRow()'s effectivePitchY() is
        // the right scale for interior row-to-row spacing, but the +/-0.5
        // pixel_tolerance margin at the sensor's outermost edge represents
        // half a single pixel's own width, not half the gap-inflated row
        // spacing. Check the true physical extent (sizeY()) directly
        // instead, with the tolerance margin scaled by the real per-pixel
        // pitch_y.
        double y_half = sizeY() / 2.0 - pixel_tolerance * pitch_y;
        if (local_y < -y_half || local_y > y_half) {
            return false;
        }
        return true;
    }

    bool DetectorGeo::hitMasked(double local_x, double local_y, int tolerance) const {
        int row = static_cast<int>(std::floor(getRow(local_y) + 0.5));
        int column = static_cast<int>(std::floor(getColumn(local_x) + 0.5));
        for (int r = row - tolerance; r <= row + tolerance; ++r) {
            for (int c = column - tolerance; c <= column + tolerance; ++c) {
                if (isMasked(c, r)) return true;
            }
        }
        return false;
    }

    bool DetectorGeo::isWithinROI(double local_x, double local_y) const {
        if (roi.empty()) return true;
        std::pair<int, int> probe{static_cast<int>(getColumn(local_x)), static_cast<int>(getRow(local_y))};
        return windingNumber(probe, roi) != 0;
    }

    void DetectorGeo::updateMatrices() {
        Eigen::AngleAxisd rotX(orientation.x(), Eigen::Vector3d::UnitX());
        Eigen::AngleAxisd rotY(orientation.y(), Eigen::Vector3d::UnitY());
        Eigen::AngleAxisd rotZ(orientation.z(), Eigen::Vector3d::UnitZ());
        R = (rotZ * rotY * rotX).matrix();
        R_inv = R.transpose();

        if (spatial_res_x > 0 && spatial_res_y > 0) {
            Eigen::Matrix3d V_loc = Eigen::Matrix3d::Zero();
            V_loc(0,0) = spatial_res_x * spatial_res_x;
            V_loc(1,1) = spatial_res_y * spatial_res_y;
            V_loc(2,2) = 1e-9;
            
            Eigen::Matrix3d V_glob = R * V_loc * R_inv;
            Eigen::Matrix2d V;
            V << V_glob(0,0), V_glob(0,1), V_glob(1,0), V_glob(1,1);
            V_inv = V.inverse();
        } else {
            V_inv = Eigen::Matrix2d::Identity();
        }
    }

    bool DetectorGeo::isPixelDetector() const {
        return n_pixels_x > 0 && n_pixels_y > 0;
    }

    GeometryManager& GeometryManager::getInstance() {
        static GeometryManager instance;
        return instance;
    }

    void GeometryManager::loadGeometry(const std::string& path) {
        toml::table tbl;
        try {
            tbl = toml::parse_file(path);
        } catch (const toml::parse_error& err) {
            throw std::runtime_error("Geometry TOML parse error in " + path + ": " + std::string(err.description()));
        }

        // Preserve file declaration order (source line) rather than
        // whatever order toml++'s own table iteration happens to give -
        // some downstream iteration/log ordering depends on it, matching
        // the old line-scanner's natural behavior.
        struct Entry {
            std::string name;
            const toml::table* table;
            uint32_t line;
        };
        std::vector<Entry> entries;
        for (auto&& [key, value] : tbl) {
            if (key == "General") continue;  // file-level metadata, not a detector - see below
            if (auto t = value.as_table()) entries.push_back({std::string(key.str()), t, value.source().begin.line});
        }
        std::sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b) { return a.line < b.line; });

        for (auto const& e : entries) {
            Configuration cfg(e.name, *e.table);
            DetectorGeo& d = detectors_[e.name];
            d.name = e.name;

            // Capture key order from *e.table (still has valid
            // source-position info at this point) before it's copied
            // into d.raw_table below - see raw_key_order's own docs for
            // why this can't be recovered later from the copy.
            {
                std::vector<std::pair<std::string, uint32_t>> keyed;
                for (auto&& [k, v] : *e.table) keyed.push_back({std::string(k.str()), v.source().begin.line});
                std::sort(keyed.begin(), keyed.end(), [](auto const& a, auto const& b) { return a.second < b.second; });
                for (auto const& [k, line] : keyed) d.raw_key_order.push_back(k);
            }
            d.raw_table = *e.table;
            detector_names_.push_back(e.name);

            d.type = cfg.get<std::string>("type", "");
            d.role = cfg.get<std::string>("role", "");
            if (d.role == "reference") { d.is_reference = true; reference_name_ = e.name; }

            auto pos = cfg.getArray<double>("position");
            if (pos.size() >= 3) d.position = Eigen::Vector3d(pos[0], pos[1], pos[2]);
            auto ori = cfg.getArray<double>("orientation");
            if (ori.size() >= 3) d.orientation = Eigen::Vector3d(ori[0], ori[1], ori[2]);

            auto pitch = cfg.getArray<double>("pixel_pitch");
            if (pitch.size() >= 2) { d.pitch_x = pitch[0]; d.pitch_y = pitch[1]; }

            d.pitch_x_gap = cfg.get<double>("pitch_x_gap", 0.0);
            d.pitch_y_gap = cfg.get<double>("pitch_y_gap", 0.0);

            auto npix = cfg.getArray<int>("number_of_pixels");
            if (npix.size() >= 2) { d.n_pixels_x = npix[0]; d.n_pixels_y = npix[1]; }

            auto sres = cfg.getArray<double>("spatial_resolution");
            if (sres.size() >= 2) { d.spatial_res_x = sres[0]; d.spatial_res_y = sres[1]; }

            d.layout = cfg.get<std::string>("layout", "");
            if (d.layout == "OrthogonalLayout") d.layout_impl = &OrthogonalLayout::instance();
            else d.layout_impl = &GenericLayout::instance();

            d.source_sensor = cfg.get<std::string>("source_sensor", "");
            d.source_channels = cfg.getArray<int>("source_channels");

            // R/R_inv/V_inv default to identity and are never derived from
            // the orientation just parsed above without this call; any
            // detector loaded with a non-zero tilt would have Track::fit()'s
            // rotated-plane path (and anything else reading det.R directly)
            // silently use identity instead of its real rotation until some
            // alignment module happened to call updateMatrices() first.
            d.updateMatrices();

            // Shared across every module (ClusteringSpatial,
            // AnalysisEfficiency, ...) instead of each one loading/parsing
            // its own copy.
            std::string mask_path = cfg.get<std::string>("mask_file", "");
            if (!mask_path.empty()) {
                std::ifstream mfile(mask_path);
                if (mfile.is_open()) {
                    std::string type;
                    int col, row;
                    while (mfile >> type >> col >> row) {
                        if (type == "p" || type == "c") d.masked_pixels.insert({col, row});
                    }
                }
            }

            // roi = [col1,row1,col2,row2,...] (flat, polygon vertices in
            // pixel coordinates) - not used by any current production
            // geometry, defaults to empty (= no restriction).
            auto roi_flat = cfg.getArray<int>("roi");
            for (size_t i = 0; i + 1 < roi_flat.size(); i += 2) d.roi.push_back({roi_flat[i], roi_flat[i + 1]});
        }

        // A geometry file self-describes where its own dies live (if
        // anywhere) - every pipeline stage that loads this file gets dies
        // for free, with no per-config opt-in required.
        if (auto general = tbl.get_as<toml::table>("General")) {
            Configuration gen_cfg("General", *general);
            dies_file_ = gen_cfg.get<std::string>("dies_file", "");
            if (!dies_file_.empty()) loadDies(dies_file_);
        }
    }

    void GeometryManager::loadDies(const std::string& path) {
        toml::table tbl;
        try {
            tbl = toml::parse_file(path);
        } catch (const toml::parse_error& err) {
            throw std::runtime_error("Dies TOML parse error in " + path + ": " + std::string(err.description()));
        }

        for (auto&& [key, value] : tbl) {
            std::string name(key.str());
            if (!hasDetector(name)) continue;
            auto t = value.as_table();
            if (!t) continue;
            Configuration cfg(name, *t);

            auto& d = detectors_.at(name);
            d.custom_dies.clear();
            for (auto& die_cfg : cfg.getArrayOfTables("dies")) {
                DetectorGeo::CustomDie die;
                die.name = die_cfg.get<std::string>("name", "");
                for (auto const& pair : die_cfg.getNestedArray<int>("pixels")) {
                    if (pair.size() >= 2) die.pixels.insert({pair[0], pair[1]});
                }
                d.custom_dies.push_back(std::move(die));
            }
        }
    }

    void GeometryManager::saveGeometry(const std::string& path) {
        std::ofstream file(path);
        if (!file.is_open()) throw std::runtime_error("Could not write to " + path);

        // toml::table's own iteration order (and hence operator<<'s
        // default formatting) is alphabetical, not insertion order -
        // verified directly against this toml++ version, true even for a
        // freshly-populated table. Detector order is reconstructed from
        // detector_names_ (already correctly ordered - see
        // loadGeometry()); each detector's own key order comes from
        // DetectorGeo::raw_key_order, captured at parse time since
        // raw_table's own copy has lost that info (see its own docs).

        // [General] first, matching the migration script's own
        // convention - order between it and the detector sections
        // doesn't matter for correctness, just readability.
        if (!dies_file_.empty()) {
            file << "[General]\n";
            file << "dies_file = \"" << dies_file_ << "\"\n\n";
        }

        for (const auto& name : detector_names_) {
            auto& d = detectors_.at(name);
            const toml::table& t = d.raw_table;

            // Best-unit (um/mm) + default-precision(6) formatting, matching
            // corry's Configuration::set(key, val, units)/Units::display
            // exactly - the only fields any alignment module actually
            // mutates; everything else round-trips verbatim from raw_table.
            toml::array position_arr{
                Units::display(d.position.x(), {"um", "mm"}),
                Units::display(d.position.y(), {"um", "mm"}),
                Units::display(d.position.z(), {"um", "mm"})
            };
            toml::array orientation_arr{
                Units::display(d.orientation.x(), {"deg"}),
                Units::display(d.orientation.y(), {"deg"}),
                Units::display(d.orientation.z(), {"deg"})
            };

            file << "[" << name << "]\n";
            for (auto const& key : d.raw_key_order) {
                file << key << " = ";
                if (key == "position" || key == "orientation") {
                    // GCC's -Warray-bounds misfires here: with the array
                    // formatter inlined into operator<<, it mistakes
                    // toml::array's internal table/array variant check for
                    // an out-of-bounds access into these locals - a known
                    // false positive in toml++, not a real one, verified by
                    // there being nothing array-indexed in this scope at all.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Warray-bounds"
#endif
                    file << (key == "position" ? position_arr : orientation_arr) << "\n";
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic pop
#endif
                } else {
                    file << toml::node_view<const toml::node>(t.get(key)) << "\n";
                }
            }
            file << "\n";
        }
    }

    bool GeometryManager::hasDetector(const std::string& name) const { return detectors_.find(name) != detectors_.end(); }
    DetectorGeo& GeometryManager::getDetector(const std::string& name) { return detectors_.at(name); }
    std::vector<std::string> GeometryManager::getDetectorNames() const { return detector_names_; }
    std::string GeometryManager::getReferenceName() const { return reference_name_; }
}