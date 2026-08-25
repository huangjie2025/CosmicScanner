# Re-trigger CMake
cmake_minimum_required(VERSION 3.15)
project(CosmicCommon CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB COMMON_SOURCES "src/*.cpp")
file(GLOB COMMON_HEADERS "src/*.h")

add_library(CosmicCommon STATIC ${COMMON_SOURCES})

target_include_directories(CosmicCommon PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Dependencies - Use local Eigen3 directly (resolve symlink for compiler compatibility)
set(EIGEN3_INCLUDE_DIR "/disk228/huangjie/PanOmics-v1.0-release/eigen3/eigen-3.4.0" CACHE PATH "Eigen3 include dir")
set(EIGEN3_FOUND TRUE CACHE BOOL "Eigen3 found")
target_include_directories(CosmicCommon PUBLIC "/disk228/huangjie/PanOmics-v1.0-release/eigen3/eigen-3.4.0")
add_library(Eigen3::Eigen INTERFACE IMPORTED)
set_target_properties(Eigen3::Eigen PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "/disk228/huangjie/PanOmics-v1.0-release/eigen3/eigen-3.4.0")
target_link_libraries(CosmicCommon PUBLIC Eigen3::Eigen)

# ZLIB and zstd from local deps
find_path(ZLIB_INCLUDE_DIR NAMES zlib.h zconf.h PATHS "/disk228/huangjie/PanOmics-v1.0-release/local_deps/include" NO_DEFAULT_PATH)
find_library(ZLIB_LIBRARY NAMES z PATHS "/disk228/huangjie/PanOmics-v1.0-release/local_deps/lib" NO_DEFAULT_PATH)
if(ZLIB_INCLUDE_DIR AND ZLIB_LIBRARY)
    set(ZLIB_FOUND TRUE)
    add_library(ZLIB::ZLIB INTERFACE IMPORTED)
    set_target_properties(ZLIB::ZLIB PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${ZLIB_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${ZLIB_LIBRARY}"
    )
    target_link_libraries(CosmicCommon PUBLIC ZLIB::ZLIB)
else()
    find_package(ZLIB REQUIRED)
    target_link_libraries(CosmicCommon PUBLIC ZLIB::ZLIB)
endif()

# zstd from local deps
find_path(ZSTD_INCLUDE_DIR NAMES zstd.h PATHS "/disk228/huangjie/PanOmics-v1.0-release/local_deps/include" NO_DEFAULT_PATH)
find_library(ZSTD_LIBRARY NAMES zstd PATHS "/disk228/huangjie/PanOmics-v1.0-release/local_deps/lib" NO_DEFAULT_PATH)
if(ZSTD_INCLUDE_DIR AND ZSTD_LIBRARY)
    set(ZSTD_FOUND TRUE)
    add_library(zstd::libzstd INTERFACE IMPORTED)
    set_target_properties(zstd::libzstd PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${ZSTD_INCLUDE_DIR}"
        INTERFACE_LINK_LIBRARIES "${ZSTD_LIBRARY}"
    )
    target_link_libraries(CosmicCommon PUBLIC zstd::libzstd)
else()
    find_package(zstd REQUIRED)
    target_link_libraries(CosmicCommon PUBLIC zstd::libzstd)
endif()

# OpenMP
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(CosmicCommon PUBLIC OpenMP::OpenMP_CXX)
endif()

# Trigger re-glob
# Added allele_harmonizer.cpp and genotype_converter.cpp
# Added spa_test.cpp and ridge_regenie.cpp
# Added rare_variant_test.cpp
# Added qc_engine.cpp
# Added slq_estimator.cpp
# Files will be picked up by GLOB