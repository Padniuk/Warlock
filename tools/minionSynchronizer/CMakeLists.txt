# --- minionSynchronizer CMake Configuration ---

add_executable(minionSynchronizer
    main.cpp
    Synchronizer.cpp
)

target_include_directories(minionSynchronizer PRIVATE
    ${CMAKE_SOURCE_DIR}/src/
    ${HDF5_INCLUDE_DIRS}
)

target_link_libraries(minionSynchronizer PRIVATE
    WarlockCore
    ${HDF5_CXX_LIBRARIES}
)

install(TARGETS minionSynchronizer RUNTIME DESTINATION bin)
