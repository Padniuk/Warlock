# --- minionRawConverter CMake Configuration ---

add_executable(minionRawConverter
    main.cpp
    RawConverter.cpp
    AppendableDataset.cpp
)

# IDEAL STATE: Set include paths for Warlock core and external dependencies
target_include_directories(minionRawConverter PRIVATE
    ${CMAKE_SOURCE_DIR}/src/
    ${EUDAQ_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
)

# IDEAL STATE: Provide hints to the linker for EUDAQ libraries
# Based on your environment, EUDAQ libs are in /eudaq/eudaq/lib
target_link_directories(minionRawConverter PRIVATE
    /eudaq/eudaq/lib
    ${EUDAQ_LIBRARY_DIRS}
)

# Link against the framework core and the external libraries
target_link_libraries(minionRawConverter PRIVATE
    WarlockCore
    ${EUDAQ_LIBRARIES}
    ${HDF5_CXX_LIBRARIES}
)

# Installation logic to the 'bin' folder
install(TARGETS minionRawConverter RUNTIME DESTINATION bin)