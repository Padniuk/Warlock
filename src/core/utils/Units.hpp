/**
 * @file Units.hpp
 * @brief Parses/formats physical quantities with unit suffixes (e.g. `"120.0GeV"`).
 */

#ifndef WARLOCK_UNITS_HPP
#define WARLOCK_UNITS_HPP

#include <string>
#include <map>
#include <vector>

namespace framework {

    /**
     * @brief Static conversion table between a config file's unit-suffixed
     * strings and Warlock's own internal units (mm for length, rad for
     * angle, ...). See Configuration::get() for where this is actually
     * invoked from config parsing.
     */
    class Units {
    public:
        /// @return `str` (e.g. `"0.05mm"`, `"120.0GeV"`) parsed and
        /// converted to Warlock's internal unit for that quantity.
        static double get(const std::string& str);
        /// @return `val` formatted with its own default display unit.
        static std::string display(double val);

        /**
         * @brief Formats `val` (an internal value in mm or rad) using
         * whichever unit in `units` displays it best, matching corry's own
         * `Units::display` algorithm exactly: prefers the unit whose
         * converted value is >= 1 and closest to 1, falling back to the
         * unit closest to 1 from below if none reach it. Used for
         * corry-parity geometry-file output.
         */
        static std::string display(double val, const std::vector<std::string>& units);
    private:
        static void initialize();
        static std::map<std::string, double> units_;
        static bool initialized_;
    };
}
#endif
