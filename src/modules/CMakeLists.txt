file(GLOB_RECURSE MODULE_SOURCES "*.cpp")

add_library(WarlockModules OBJECT ${MODULE_SOURCES})

target_include_directories(WarlockModules PUBLIC ${PROJECT_SOURCE_DIR}/src)
target_link_libraries(WarlockModules PUBLIC 
    WarlockObjects 
    Eigen3::Eigen 
    tomlplusplus::tomlplusplus 
    ${HDF5_CXX_LIBRARIES}
)

# Pass the object library to the main scope
set(WARLOCK_MODULE_LIBRARIES WarlockModules PARENT_SCOPE)