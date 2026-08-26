#include "core/Module.hpp"

namespace framework {
    Module::Module(Configuration cfg, Configuration global_cfg, ThreadPool* pool)
        : config(std::move(cfg)), global_config(std::move(global_cfg)), thread_pool(pool) {}

    std::string Module::getName() const {
        return config.getName();
    }

    void Module::setOutputPath(std::filesystem::path path) {
        output_path_ = std::move(path);
    }

    std::filesystem::path Module::getOutputPath() const {
        return output_path_;
    }
}