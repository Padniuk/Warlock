# Macro to define the build option for a module
macro(WARLOCK_MODULE name)
    option(BUILD_${name} "Build module ${name}?" ON)
endmacro()

# Macro to create a library for the module and register it
macro(WARLOCK_MODULE_SOURCES name)
    if(BUILD_${name})
        # Create a static library for this specific module
        add_library(${name} STATIC ${ARGN})
        
        # All modules need access to the src/ directory and dependencies
        target_include_directories(${name} PUBLIC ${PROJECT_SOURCE_DIR}/src)
        target_link_libraries(${name} PUBLIC 
            WarlockObjects 
            Eigen3::Eigen 
            tomlplusplus::tomlplusplus 
            ${HDF5_CXX_LIBRARIES}
        )

        # Append this library to a global list so the main executable can link it
        set(WARLOCK_MODULE_LIBRARIES ${WARLOCK_MODULE_LIBRARIES} ${name} CACHE INTERNAL "Module libraries")
    endif()
endmacro()

# Macro to handle installation of the module
macro(WARLOCK_MODULE_INSTALL name)
    if(BUILD_${name})
        install(TARGETS ${name} DESTINATION lib)
    endif()
endmacro()