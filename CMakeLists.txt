cmake_minimum_required(VERSION 3.14)
project(Warlock VERSION 1.0.0)

# --- 1. Basic Setup ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/install" CACHE PATH "..." FORCE)
endif()

include(FetchContent)
include(GNUInstallDirs)

# M_PI etc. are POSIX/GNU extensions that GCC/Clang expose in <cmath> by
# default, but MSVC only defines them when this macro is set before the
# header is first included - a command-line -D applies before any file is
# processed, so it works regardless of header include order. Harmless no-op
# on GCC/Clang.
add_compile_definitions(_USE_MATH_DEFINES)

# --- 2. Build Options ---
# ROOT is notoriously hard to get working on Windows and isn't needed by the
# core `warlock` executable or any of its modules - only by this optional
# plotting tool. Auto-detect it so a plain `cmake ..` on a machine without
# ROOT (e.g. a fresh Windows setup) configures cleanly with the tool simply
# skipped, rather than hard-failing; explicitly passing
# -DBUILD_MINION_PLOTTER=ON still fails loudly if ROOT truly isn't there.
find_package(ROOT QUIET COMPONENTS Hist RIO Core)
option(BUILD_MINION_PLOTTER   "Build HDF5-to-ROOT tool? (requires ROOT)" ${ROOT_FOUND})
option(BUILD_MINION_CONVERTER "Build EUDAQ-to-HDF5 tool?" ON)
option(BUILD_MINION_SYNCHRONIZER "Build CAEN/MIMOSA desync detector+fixer tool? (needs only HDF5, already required)" ON)

# --- 3. Dependencies ---
FetchContent_Declare(tomlplusplus GIT_REPOSITORY https://github.com/marzer/tomlplusplus.git GIT_TAG v3.3.0)
FetchContent_MakeAvailable(tomlplusplus)

# Eigen is header-only for everything Warlock actually uses (dense matrix/
# vector ops) - but its own CMakeLists.txt is written for building EIGEN
# ITSELF as a full project (unit tests, docs, benchmarks, sparse-solver
# backend detection: CHOLMOD/UMFPACK/KLU/SuperLU/PaStiX/SPQR/MPFR/FFTW/
# GoogleHash/Adolc, CUDA/Fortran/Qt4 probing...). FetchContent_MakeAvailable
# would add_subdirectory() straight into that, which both (a) floods
# ccmake/cmake-gui with hundreds of EIGEN_*/unrelated options nothing here
# ever touches, and (b) inflates `make`'s file-count estimate with Eigen's
# own internal targets, which is why the progress percentage used to jump in
# a few big, uneven steps instead of climbing smoothly through Warlock's own
# ~30 files. FetchContent_Populate() only fetches the source - Eigen's
# CMakeLists.txt is never processed - and a plain INTERFACE target aliased
# to the name everything already links against (Eigen3::Eigen) is all a
# header-only dependency needs.
FetchContent_Declare(Eigen3 GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git GIT_TAG 3.4.0)
FetchContent_GetProperties(Eigen3)
if(NOT eigen3_POPULATED)
    FetchContent_Populate(Eigen3)
endif()
add_library(eigen_headers INTERFACE)
target_include_directories(eigen_headers INTERFACE ${eigen3_SOURCE_DIR})
add_library(Eigen3::Eigen ALIAS eigen_headers)

find_package(HDF5 COMPONENTS CXX REQUIRED)
find_package(Threads REQUIRED)

if(BUILD_MINION_PLOTTER AND NOT ROOT_FOUND)
    find_package(ROOT REQUIRED COMPONENTS Hist RIO Core)
endif()

if(BUILD_MINION_CONVERTER)
    # Hint for find_package if needed
    set(eudaq_DIR "/eudaq/eudaq/cmake" CACHE PATH "Path to EUDAQ cmake files")
    find_package(eudaq REQUIRED)
    message(STATUS "Found EUDAQ: ${EUDAQ_INCLUDE_DIRS}")
endif()

# All of the above is real, required configuration - HDF5, and ROOT when
# BUILD_MINION_PLOTTER needs it, aren't optional; the software won't build
# without them. mark_as_advanced() doesn't remove or hide that fact, it only
# keeps FetchContent's own bookkeeping and each dependency's individually-
# discovered library paths out of ccmake/cmake-gui's DEFAULT view, which is
# the only thing this changes - all of it is still fully visible (and
# editable, if some path genuinely needs overriding) in advanced mode
# (ccmake: press 't'; cmake-gui: the "Advanced" checkbox; command line:
# `cmake -LA`). A plain `ccmake ..` is then just BUILD_MINION_CONVERTER,
# BUILD_MINION_PLOTTER, and CMake's own two standard entries
# (CMAKE_BUILD_TYPE, CMAKE_INSTALL_PREFIX) - the knobs actually meant to be
# turned, without pretending the rest doesn't exist.
mark_as_advanced(
    FETCHCONTENT_BASE_DIR
    FETCHCONTENT_FULLY_DISCONNECTED
    FETCHCONTENT_QUIET
    FETCHCONTENT_SOURCE_DIR_EIGEN3
    FETCHCONTENT_SOURCE_DIR_TOMLPLUSPLUS
    FETCHCONTENT_UPDATES_DISCONNECTED
    FETCHCONTENT_UPDATES_DISCONNECTED_EIGEN3
    FETCHCONTENT_UPDATES_DISCONNECTED_TOMLPLUSPLUS
    HDF5_CXX_LIBRARY_crypto
    HDF5_CXX_LIBRARY_curl
    HDF5_CXX_LIBRARY_dl
    HDF5_CXX_LIBRARY_hdf5
    HDF5_CXX_LIBRARY_hdf5_cpp
    HDF5_CXX_LIBRARY_m
    HDF5_CXX_LIBRARY_pthread
    HDF5_CXX_LIBRARY_sz
    HDF5_CXX_LIBRARY_z
    ROOT_DIR
    ROOT_genmap_CMD
    ROOT_rootdraw_CMD
)

# --- 4. Core Targets ---

# Object headers interface
add_library(WarlockObjects INTERFACE)
target_include_directories(WarlockObjects INTERFACE src/)

# Main Framework Engine
file(GLOB_RECURSE CORE_SOURCES "src/core/*.cpp")
add_library(WarlockCore STATIC ${CORE_SOURCES})
target_link_libraries(WarlockCore PUBLIC
    ${HDF5_CXX_LIBRARIES}
    Threads::Threads
    tomlplusplus::tomlplusplus
    Eigen3::Eigen
    WarlockObjects
)

# --- 5. Modules ---
# Initialize the list of libraries to link
set(WARLOCK_MODULE_LIBRARIES "" CACHE INTERNAL "")
add_subdirectory(src/modules)

# --- 6. Main Executable ---
add_executable(warlock main.cpp)
target_link_libraries(warlock PRIVATE WarlockCore ${WARLOCK_MODULE_LIBRARIES})

# Compiler Flags
# Note: deliberately no -march=native - it bakes in instructions specific to
# the exact CPU that built the binary, so a binary built on one machine can
# crash with an illegal-instruction fault on another (even same OS/arch), and
# it isn't available on MSVC at all. Not worth the deployment portability risk.
if(NOT MSVC)
    target_compile_options(WarlockCore PRIVATE -Wall -Wextra -O3)
    target_compile_options(warlock PRIVATE -Wall -Wextra -O3)
    # WarlockModules (every module - WaveformProcessingCAEN, Tracking4D,
    # ClusteringSpatial, EventLoaderHDF5, ...) never had this applied at
    # all: its own flags.make showed only "-std=gnu++17", no optimization
    # flag whatsoever, meaning every module has been compiled at the
    # compiler's unoptimized default this entire time - no inlining, no
    # vectorization of the Eigen expressions used throughout the per-hit/
    # per-waveform/per-cluster hot loops, nothing. Found while investigating
    # whether manually marking specific functions `inline` was worth doing;
    # it wasn't the real issue - the module target just wasn't being
    # optimized at all.
    target_compile_options(WarlockModules PRIVATE -Wall -Wextra -O3)
else()
    target_compile_options(WarlockCore PRIVATE /O2)
    target_compile_options(warlock PRIVATE /O2)
    target_compile_options(WarlockModules PRIVATE /O2)
endif()

# --- 7. Tools & Installation ---
if(BUILD_MINION_PLOTTER)
    add_subdirectory(tools/minionRootPlotter)
endif()
if(BUILD_MINION_CONVERTER)
    add_subdirectory(tools/minionRawConverter)
endif()
if(BUILD_MINION_SYNCHRONIZER)
    add_subdirectory(tools/minionSynchronizer)
endif()

install(TARGETS warlock RUNTIME DESTINATION bin)