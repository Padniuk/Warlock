#include "core/config/Configuration.hpp"

namespace framework {
    Configuration::Configuration(std::string name, toml::table table)
        : name_(std::move(name)), table_(std::move(table)) {}
}