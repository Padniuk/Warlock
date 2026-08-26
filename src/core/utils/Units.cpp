#include "core/utils/Units.hpp"
#include <regex>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <limits>

namespace framework {
    std::map<std::string, double> Units::units_;
    bool Units::initialized_ = false;

    void Units::initialize() {
        if (initialized_) return;
        units_["nm"] = 1e-6; units_["um"] = 1e-3; units_["mm"] = 1.0; units_["cm"] = 10.0; units_["m"] = 1000.0;
        units_["ps"] = 1e-3; units_["ns"] = 1.0; units_["us"] = 1e3; units_["ms"] = 1e6; units_["s"] = 1e9;
        units_["deg"] = M_PI / 180.0; units_["rad"] = 1.0; units_["mrad"] = 0.001;
        initialized_ = true;
    }

    double Units::get(const std::string& str) {
        initialize();
        std::regex re("([-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?)\\s*([a-zA-Z]+)?");
        std::smatch match;
        if (std::regex_search(str, match, re)) {
            double val = std::stod(match[1]);
            std::string unit = match[2];
            if (unit.empty()) return val;
            if (units_.count(unit)) return val * units_[unit];
        }
        return 0.0;
    }

    std::string Units::display(double val) {
        std::ostringstream oss;
        if (std::abs(val) < 1.0 && val != 0) oss << std::fixed << std::setprecision(2) << val * 1000.0 << "um";
        else oss << std::fixed << std::setprecision(4) << val << "mm";
        return oss.str();
    }

    std::string Units::display(double val, const std::vector<std::string>& units) {
        initialize();
        if (units.empty()) return display(val);

        int best_exponent = std::numeric_limits<int>::min();
        std::string best_unit = units.front();
        for (const auto& unit : units) {
            auto it = units_.find(unit);
            if (it == units_.end() || it->second == 0.0) continue;
            double value = val / it->second;
            int exponent = 0;
            std::frexp(value, &exponent);
            if ((best_exponent <= 0 && exponent > best_exponent) || (exponent > 0 && exponent < best_exponent)) {
                best_exponent = exponent;
                best_unit = unit;
            }
        }

        std::ostringstream oss;
        oss << (val / units_[best_unit]) << best_unit;
        return oss.str();
    }
}