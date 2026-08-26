add_executable(minionRootPlotter
    main.cpp
    H5ToRootConverter.cpp
)

target_include_directories(minionRootPlotter PRIVATE
    ${CMAKE_SOURCE_DIR}/src/
    ${ROOT_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
)

target_link_libraries(minionRootPlotter PRIVATE
    ${HDF5_CXX_LIBRARIES}
    ROOT::Hist
    ROOT::RIO
    ROOT::Core
)