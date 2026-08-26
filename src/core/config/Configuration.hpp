/**
 * @file Configuration.hpp
 * @brief Typed accessor for a single module's TOML configuration block.
 */

#ifndef WARLOCK_CONFIGURATION_HPP
#define WARLOCK_CONFIGURATION_HPP

#include <string>
#include <vector>
#include <optional>
#include <type_traits>
#include <toml++/toml.hpp>
#include "core/utils/Units.hpp"

namespace framework {

    /**
     * @brief Read-only, typed view over one TOML table (either the
     * framework-wide `[Warlock]` block, or a single module's own block).
     */
    class Configuration {
    public:
        Configuration() = default;
        /// @param name The block/module instance name this config belongs to.
        /// @param table The already-parsed TOML table for this block.
        Configuration(std::string name, toml::table table);

        /**
         * @brief Reads `key` as type `T`, or `default_value` if absent.
         *
         * For arithmetic `T` (other than `bool`), a value given as a
         * unit-suffixed string (e.g. `"120.0GeV"`, `"0.05mm"`) is parsed
         * via Units::get() instead of failing - lets config files write
         * physical quantities in whatever unit is natural, not just plain
         * numbers.
         */
        template <typename T>
        T get(const std::string& key, T default_value) const {
            auto node = table_.get(key);
            if (!node) return default_value;
            if (auto val = node->value<T>()) return *val;
            if constexpr (std::is_arithmetic_v<T> && !std::is_same_v<T, bool>) {
                if (auto str_val = node->value<std::string>()) {
                    return static_cast<T>(Units::get(*str_val));
                }
            }
            return default_value;
        }

        /**
         * @brief Reads `key` as an array of `T`.
         *
         * Accepts either a genuine TOML array, or a single bare value
         * (treated as a one-element array) - lets a config write
         * `cfd = 0.7` instead of `cfd = [0.7]` when only one value is
         * needed. Missing keys yield an empty vector. Same
         * unit-suffixed-string fallback as get() applies per element.
         */
        template <typename T>
        std::vector<T> getArray(const std::string& key) const {
            std::vector<T> result;
            auto node = table_.get(key);
            if (!node) return result;

            if (node->is_array()) {
                auto& arr = *node->as_array();
                for (auto&& el : arr) {
                    if (auto val = el.value<T>()) result.push_back(*val);
                    else if constexpr (std::is_arithmetic_v<T> && !std::is_same_v<T, bool>) {
                        if (auto str_val = el.value<std::string>()) result.push_back(static_cast<T>(Units::get(*str_val)));
                    }
                }
            } else {
                if (auto val = node->value<T>()) result.push_back(*val);
                else if constexpr (std::is_arithmetic_v<T> && !std::is_same_v<T, bool>) {
                    if (auto str_val = node->value<std::string>()) result.push_back(static_cast<T>(Units::get(*str_val)));
                }
            }
            return result;
        }

        /**
         * @brief Reads `key` as a genuine array of arrays (e.g.
         * `pixels = [[0,0], [1,0]]`), each inner array read with the same
         * per-element logic as getArray<T>().
         *
         * Missing keys, or a key that isn't an array of arrays, yield an
         * empty vector. Unlike getArray<T>(), there's no "bare scalar"
         * fallback here - an inner element that isn't a real array is
         * skipped rather than treated as a one-element one, since a
         * mixed flat/nested list almost always means a config mistake
         * worth surfacing as a shorter-than-expected result rather than
         * silently reinterpreting it.
         */
        template <typename T>
        std::vector<std::vector<T>> getNestedArray(const std::string& key) const {
            std::vector<std::vector<T>> result;
            auto node = table_.get(key);
            if (!node || !node->is_array()) return result;

            for (auto&& el : *node->as_array()) {
                auto inner_arr = el.as_array();
                if (!inner_arr) continue;
                std::vector<T> inner;
                for (auto&& inner_el : *inner_arr) {
                    if (auto val = inner_el.value<T>()) inner.push_back(*val);
                    else if constexpr (std::is_arithmetic_v<T> && !std::is_same_v<T, bool>) {
                        if (auto str_val = inner_el.value<std::string>()) inner.push_back(static_cast<T>(Units::get(*str_val)));
                    }
                }
                result.push_back(std::move(inner));
            }
            return result;
        }

        /**
         * @brief Reads `key` as a list of nested tables, each wrapped as
         * its own Configuration.
         *
         * Accepts either a genuine TOML array-of-tables (`[[key]]` blocks)
         * or a single bare table (treated as a one-element array, same
         * "bare scalar = one-element array" convention as getArray()).
         * Missing keys, or a key that isn't a table/array-of-tables,
         * yield an empty vector. Lets a config express real nested
         * structure (e.g. a detector's own named dies) while reading each
         * element with the same get<T>()/getArray<T>() API as any other
         * Configuration.
         */
        std::vector<Configuration> getArrayOfTables(const std::string& key) const {
            std::vector<Configuration> result;
            auto node = table_.get(key);
            if (!node) return result;

            if (auto arr = node->as_array()) {
                for (auto&& el : *arr) {
                    if (auto t = el.as_table()) result.emplace_back(key, *t);
                }
            } else if (auto t = node->as_table()) {
                result.emplace_back(key, *t);
            }
            return result;
        }

        /// @return Whether `key` is present in this block at all.
        bool has(const std::string& key) const { return table_.contains(key); }
        /// @return The block/module instance name this config belongs to.
        std::string getName() const { return name_; }

    private:
        std::string name_;
        toml::table table_;
    };
}
#endif
