cmake_minimum_required(VERSION 3.15)
project(CosmicScanner VERSION 1.1.3 LANGUAGES CXX)

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

set(COSMICSCANNER_GIT_COMMIT_OVERRIDE "" CACHE STRING "Override Git commit embedded in scan metadata")
if(COSMICSCANNER_GIT_COMMIT_OVERRIDE)
  set(COSMICSCANNER_GIT_COMMIT "${COSMICSCANNER_GIT_COMMIT_OVERRIDE}")
elseif(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/COSMICSCANNER_GIT_COMMIT")
  file(READ "${CMAKE_CURRENT_SOURCE_DIR}/COSMICSCANNER_GIT_COMMIT" COSMICSCANNER_GIT_COMMIT)
  string(STRIP "${COSMICSCANNER_GIT_COMMIT}" COSMICSCANNER_GIT_COMMIT)
else()
  execute_process(
    COMMAND git rev-parse --short=12 HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE COSMICSCANNER_GIT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
  )
endif()
if(NOT COSMICSCANNER_GIT_COMMIT)
  set(COSMICSCANNER_GIT_COMMIT "unknown")
endif()

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
if(WIN32)
  target_link_libraries(cosmicscanner PRIVATE psapi)
endif()
target_compile_definitions(cosmicscanner PRIVATE ${COSMICSCANNER_EXTRA_DEFS})
target_compile_definitions(cosmicscanner PRIVATE COSMICSCANNER_VERSION="${PROJECT_VERSION}")
target_compile_definitions(cosmicscanner PRIVATE COSMICSCANNER_GIT_COMMIT="${COSMICSCANNER_GIT_COMMIT}")
if(APPLE)
  set_target_properties(cosmicscanner PROPERTIES INSTALL_RPATH "@loader_path/../lib")
endif()

install(TARGETS cosmicscanner RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
install(FILES README.md LICENSE CHANGELOG.md CITATION.cff DESTINATION ${CMAKE_INSTALL_DATADIR}/cosmicscanner)
install(DIRECTORY docs/ DESTINATION ${CMAKE_INSTALL_DATADIR}/cosmicscanner/docs)

include(CTest)
if(BUILD_TESTING)
  add_executable(test_mv_em_update ${CMAKE_CURRENT_SOURCE_DIR}/test/test_mv_em_update.cpp)
  target_include_directories(test_mv_em_update PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries(test_mv_em_update PRIVATE Eigen3::Eigen)
  add_test(NAME mv_em_update_scaling COMMAND test_mv_em_update)

  add_executable(test_mv_reml_dense ${CMAKE_CURRENT_SOURCE_DIR}/test/test_mv_reml_dense.cpp)
  target_include_directories(test_mv_reml_dense PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  target_link_libraries(test_mv_reml_dense PRIVATE Eigen3::Eigen)
  add_test(NAME mv_reml_objective_dense COMMAND test_mv_reml_dense objective)
  add_test(NAME mv_reml_gradient_finite_difference COMMAND test_mv_reml_dense gradient)
  add_test(NAME mv_reml_update_monotonicity COMMAND test_mv_reml_dense monotonicity)

  find_program(BASH_EXECUTABLE NAMES bash)
  if(BASH_EXECUTABLE)
    add_test(
      NAME critical_regressions
      COMMAND ${BASH_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test/run_critical_regressions.sh $<TARGET_FILE:cosmicscanner>
    )
  else()
    message(WARNING "Bash not found; critical regression tests are unavailable")
  endif()
endif()
