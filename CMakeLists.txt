cmake_minimum_required(VERSION 3.15)
project(CosmicScanner VERSION 1.1.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(GNUInstallDirs)

find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)
find_package(ZLIB REQUIRED)
find_package(zstd REQUIRED)
if(TARGET zstd::libzstd)
  set(COSMICSCANNER_ZSTD_TARGET zstd::libzstd)
elseif(TARGET zstd::libzstd_shared)
  set(COSMICSCANNER_ZSTD_TARGET zstd::libzstd_shared)
elseif(TARGET zstd::libzstd_static)
  set(COSMICSCANNER_ZSTD_TARGET zstd::libzstd_static)
else()
  message(FATAL_ERROR "A supported zstd CMake target was not found")
endif()
option(COSMICSCANNER_ENABLE_BLAS "Enable BLAS acceleration" OFF)

set(GWAS_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/src/a_matrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/adaptive_trace_estimator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/bgen_reader.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/covar_reader.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/dense_solver.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/design.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/eigen_vce.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/grm_builder.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gwas_output.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/gwas_scanner.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/lmm_engine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/matrix_io.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mme_builder.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/multivariate_gemma_refit.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/multivariate_model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/multivariate_output.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/mvlmm_engine.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/pgen_reader.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/pheno_reader.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/plink_reader.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/qc_utils.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/repeatability_cache.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/repeatability_model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/repeatability_output.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/repeatability_refit_validation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/repeatability_scan_exact.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/repeatability_scan_fast.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/repeatability_scan_stcg.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/repeatability_utils.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/single_step_model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/slq_estimator.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/stcg_vce.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/univariate_model.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/v_matrix.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/src/vce.cpp
)

if(COSMICSCANNER_ENABLE_BLAS)
  find_package(BLAS)
  if(BLAS_FOUND)
    list(APPEND COSMICSCANNER_EXTRA_LIBS ${BLAS_LIBRARIES})
    list(APPEND COSMICSCANNER_EXTRA_DEFS EIGEN_USE_BLAS)
  endif()
endif()

add_executable(cosmicscanner ${GWAS_SOURCES})
target_include_directories(cosmicscanner PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(cosmicscanner PRIVATE Eigen3::Eigen OpenMP::OpenMP_CXX ZLIB::ZLIB ${COSMICSCANNER_ZSTD_TARGET} ${COSMICSCANNER_EXTRA_LIBS})
target_compile_definitions(cosmicscanner PRIVATE ${COSMICSCANNER_EXTRA_DEFS})
target_compile_definitions(cosmicscanner PRIVATE COSMICSCANNER_VERSION="${PROJECT_VERSION}")

install(TARGETS cosmicscanner RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

include(CTest)
if(BUILD_TESTING)
  add_test(
    NAME critical_regressions
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/test/run_critical_regressions.sh $<TARGET_FILE:cosmicscanner>
  )
endif()
