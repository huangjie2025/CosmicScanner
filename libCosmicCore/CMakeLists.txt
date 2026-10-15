# Re-trigger CMake
cmake_minimum_required(VERSION 3.15)
project(CosmicCommon CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB COMMON_SOURCES "src/*.cpp")
file(GLOB COMMON_HEADERS "src/*.h")

add_library(CosmicCommon STATIC ${COMMON_SOURCES})

target_include_directories(CosmicCommon PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src)

# Dependencies
# Eigen3 is header-only, find it via include directory
# Use CONDA_PREFIX if set and contains eigen3, otherwise try miniconda3 base
if (DEFINED ENV{CONDA_PREFIX})
    # Check if the conda prefix has eigen3, otherwise fall back to miniconda3 base
    if (EXISTS "$ENV{CONDA_PREFIX}/include/eigen3/Eigen/Core")
        set(EIGEN3_INCLUDE_DIR "$ENV{CONDA_PREFIX}/include/eigen3")
    else()
        find_path(EIGEN3_INCLUDE_DIR Eigen/Dense PATHS /usr/include/eigen3 /usr/local/include/eigen3)
    endif()
else()
    find_path(EIGEN3_INCLUDE_DIR Eigen/Dense PATHS /usr/include/eigen3 /usr/local/include/eigen3)
endif()

if (EIGEN3_INCLUDE_DIR AND EXISTS ${EIGEN3_INCLUDE_DIR}/Eigen/Core)
    message(STATUS "Eigen3 found at: ${EIGEN3_INCLUDE_DIR}")
    target_include_directories(CosmicCommon PUBLIC ${EIGEN3_INCLUDE_DIR})
else()
    message(FATAL_ERROR "Eigen3 not found. Please install Eigen3 or set CONDA_PREFIX.")
endif()

# BLAS/LAPACK support for Eigen acceleration
find_package(BLAS)
if(BLAS_FOUND)
    target_link_libraries(CosmicCommon PUBLIC ${BLAS_LIBRARIES})
    # Use EIGEN_USE_BLAS always; cmake will find MKL, OpenBLAS, or generic BLAS
    target_compile_definitions(CosmicCommon PUBLIC EIGEN_USE_BLAS)
    message(STATUS "BLAS found: ${BLAS_LIBRARIES}")
endif()

find_package(LAPACK)
if(LAPACK_FOUND)
    target_link_libraries(CosmicCommon PUBLIC ${LAPACK_LIBRARIES})
    target_compile_definitions(CosmicCommon PUBLIC EIGEN_USE_LAPACKE)
    message(STATUS "LAPACK found: ${LAPACK_LIBRARIES}")
endif()

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(CosmicCommon PUBLIC OpenMP::OpenMP_CXX)
endif()

find_package(ZLIB REQUIRED)
target_link_libraries(CosmicCommon PUBLIC ZLIB::ZLIB)

find_package(zstd REQUIRED)
# Upstream zstd exports static/shared targets; some distributions provide
# the historical unsuffixed target instead.
if(TARGET zstd::libzstd)
    target_link_libraries(CosmicCommon PUBLIC zstd::libzstd)
elseif(TARGET zstd::libzstd_static)
    target_link_libraries(CosmicCommon PUBLIC zstd::libzstd_static)
elseif(TARGET zstd::libzstd_shared)
    target_link_libraries(CosmicCommon PUBLIC zstd::libzstd_shared)
else()
    message(FATAL_ERROR "zstd was found but exports no supported library target")
endif()

# Trigger re-glob
# Added allele_harmonizer.cpp and genotype_converter.cpp
# Added spa_test.cpp and ridge_regenie.cpp
# Added rare_variant_test.cpp
# Added qc_engine.cpp
# Added slq_estimator.cpp
# Files will be picked up by GLOB

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    include(CTest)
    if(BUILD_TESTING)
        add_executable(test_solver_planner tests/test_solver_planner.cpp)
        target_link_libraries(test_solver_planner PRIVATE CosmicCommon)
        add_test(NAME solver_planner_unit COMMAND test_solver_planner)
        set_tests_properties(solver_planner_unit PROPERTIES LABELS "unit;solver_planner" TIMEOUT 30)

        add_executable(test_analysis_dataset tests/test_analysis_dataset.cpp)
        target_link_libraries(test_analysis_dataset PRIVATE CosmicCommon)
        add_test(NAME analysis_dataset_unit COMMAND test_analysis_dataset)
        set_tests_properties(analysis_dataset_unit PROPERTIES LABELS "unit;analysis_dataset" TIMEOUT 30)

        add_executable(test_validation tests/test_validation.cpp)
        target_link_libraries(test_validation PRIVATE CosmicCommon)
        add_test(NAME validation_unit COMMAND test_validation)
        set_tests_properties(validation_unit PROPERTIES LABELS "unit;validation" TIMEOUT 30)

        add_executable(test_run_manifest tests/test_run_manifest.cpp)
        target_link_libraries(test_run_manifest PRIVATE CosmicCommon)
        add_test(NAME run_manifest_unit COMMAND test_run_manifest)
        set_tests_properties(run_manifest_unit PROPERTIES LABELS "unit;manifest" TIMEOUT 30)

        add_executable(test_reader_compatibility tests/test_reader_compatibility.cpp)
        target_link_libraries(test_reader_compatibility PRIVATE CosmicCommon)
        add_test(NAME reader_compatibility_unit COMMAND test_reader_compatibility)
        set_tests_properties(reader_compatibility_unit PROPERTIES LABELS "unit;reader" TIMEOUT 60)

        add_executable(test_bgen_layout2_reader tests/test_bgen_layout2_reader.cpp src/bgen_reader.cpp)
        target_include_directories(test_bgen_layout2_reader PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${EIGEN3_INCLUDE_DIR})
        target_link_libraries(test_bgen_layout2_reader PRIVATE ZLIB::ZLIB zstd::libzstd)
        if(OpenMP_CXX_FOUND)
            target_link_libraries(test_bgen_layout2_reader PRIVATE OpenMP::OpenMP_CXX)
        endif()
        add_test(NAME bgen_layout2_reader_unit COMMAND test_bgen_layout2_reader)
        set_tests_properties(bgen_layout2_reader_unit PROPERTIES LABELS "unit;reader;bgen" TIMEOUT 30)

        add_executable(test_eigen_vce_ad tests/test_eigen_vce_ad.cpp)
        target_link_libraries(test_eigen_vce_ad PRIVATE CosmicCommon)
        add_test(NAME eigen_vce_ad_unit COMMAND test_eigen_vce_ad)
        set_tests_properties(eigen_vce_ad_unit PROPERTIES LABELS "unit;vce;ad" TIMEOUT 120)

        add_executable(test_ad_variance_recovery tests/test_ad_variance_recovery.cpp)
        target_link_libraries(test_ad_variance_recovery PRIVATE CosmicCommon)
        add_test(NAME ad_variance_recovery_unit COMMAND test_ad_variance_recovery)
        set_tests_properties(ad_variance_recovery_unit PROPERTIES LABELS "unit;vce;ad;validation" TIMEOUT 3600)

        add_executable(test_ad_production_kd_recovery tests/test_ad_production_kd_recovery.cpp)
        target_link_libraries(test_ad_production_kd_recovery PRIVATE CosmicCommon)
        add_test(NAME ad_production_kd_recovery COMMAND test_ad_production_kd_recovery)
        set_tests_properties(ad_production_kd_recovery PROPERTIES LABELS "validation;vce;ad;grm" TIMEOUT 1800)

        add_executable(test_ad_reml_df_regression tests/test_ad_reml_df_regression.cpp)
        target_link_libraries(test_ad_reml_df_regression PRIVATE CosmicCommon)
        add_test(NAME ad_reml_df_regression COMMAND test_ad_reml_df_regression)
        set_tests_properties(ad_reml_df_regression PROPERTIES LABELS "unit;vce;ad;validation" TIMEOUT 120)

        add_executable(test_eigen_vce_d tests/test_eigen_vce_d.cpp)
        target_link_libraries(test_eigen_vce_d PRIVATE CosmicCommon)
        add_test(NAME eigen_vce_d_unit COMMAND test_eigen_vce_d)
        set_tests_properties(eigen_vce_d_unit PROPERTIES LABELS "unit;vce;dominance" TIMEOUT 120)
    endif()
endif()
