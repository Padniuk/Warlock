/**
 * @file ModuleFactory.hpp
 * @brief Runtime registry mapping a module's TOML type name to its constructor.
 */

#ifndef WARLOCK_MODULEFACTORY_HPP
#define WARLOCK_MODULEFACTORY_HPP

#include "core/Module.hpp"
#include <functional>
#include <map>
#include <memory>
#include <string>

namespace framework {

    /**
     * @brief Singleton registry of module type name -> constructor function.
     *
     * Modules register themselves via the REGISTER_MODULE macro at static
     * initialization time; FrameworkManager looks a module type up here by
     * the name used in its TOML `[TypeName]` block to construct an instance.
     */
    class ModuleFactory {
    public:
        using CreatorFunc = std::function<std::unique_ptr<Module>(Configuration, Configuration, ThreadPool*)>;

        static ModuleFactory& getInstance() {
            static ModuleFactory instance;
            return instance;
        }

        /// Registers a module type's constructor under its TOML type name.
        void registerModule(const std::string& name, CreatorFunc func) {
            registry_[name] = std::move(func);
        }

        /**
         * @brief Constructs a registered module instance by type name.
         * @throws std::runtime_error if `name` was never registered.
         */
        std::unique_ptr<Module> createModule(const std::string& name, Configuration cfg, Configuration global_cfg, ThreadPool* pool) {
            if (registry_.find(name) == registry_.end()) {
                throw std::runtime_error("Module '" + name + "' is not registered.");
            }
            return registry_[name](std::move(cfg), std::move(global_cfg), pool);
        }

    private:
        ModuleFactory() = default;
        std::map<std::string, CreatorFunc> registry_;
    };

    /**
     * @brief Registers module type `T` with ModuleFactory as a side effect
     * of constructing a static instance of this class (see REGISTER_MODULE).
     */
    template <typename T>
    class ModuleRegistrar {
    public:
        ModuleRegistrar(const std::string& name) {
            ModuleFactory::getInstance().registerModule(name, [](Configuration cfg, Configuration g_cfg, ThreadPool* p) {
                return std::make_unique<T>(std::move(cfg), std::move(g_cfg), p);
            });
        }
    };
}

/**
 * @brief Registers a Module subclass with ModuleFactory under its own class
 * name, so it can be constructed by that name from a TOML `[ClassName]`
 * block. Place once at namespace scope in the module's own .cpp file.
 */
#define REGISTER_MODULE(ClassName) \
    static framework::ModuleRegistrar<ClassName> registrar_##ClassName(#ClassName);
#endif
