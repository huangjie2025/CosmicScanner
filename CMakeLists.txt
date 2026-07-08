cmake_minimum_required(VERSION 3.15)
project(CosmicScanner VERSION 1.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include(GNUInstallDirs)
option(COSMICSCANNER_BUILD_TESTS "Build CosmicScanner test targets" ON)

# Add submodule library
add_subdirectory(libCosmicCore)

# Find necessary libraries (e.g. OpenMP)
find_package(OpenMP REQUIRED)

# --- Source Files ---
file(GLOB COMMON_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
file(GLOB MODEL_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/models/*.cpp")
file(GLOB INPUT_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/input/*.cpp")
list(APPEND COMMON_SOURCES ${MODEL_SOURCES})
list(APPEND COMMON_SOURCES ${INPUT_SOURCES})
# 过滤掉带有 main 函数的测试或入口文件，防止链接冲突
list(REMOVE_ITEM COMMON_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
list(REMOVE_ITEM COMMON_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/gwas_app.cpp")

add_executable(cosmicscanner src/main.cpp src/gwas_app.cpp ${COMMON_SOURCES})

target_include_directories(cosmicscanner PRIVATE 
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/models
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)

target_link_libraries(cosmicscanner PRIVATE CosmicCommon OpenMP::OpenMP_CXX)
target_compile_definitions(cosmicscanner PRIVATE COSMICSCANNER_VERSION="${PROJECT_VERSION}")

if(COSMICSCANNER_BUILD_TESTS)
enable_testing()

add_test(
    NAME gemma_emmvgwas_compat_smoke
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_gemma_emmvgwas_compat_smoke.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
set_tests_properties(gemma_emmvgwas_compat_smoke PROPERTIES
    ENVIRONMENT "COSMICSCANNER_ROOT=${CMAKE_CURRENT_SOURCE_DIR};COSMICSCANNER_BIN=$<TARGET_FILE:cosmicscanner>;COSMIC_SMOKE_OUT=${CMAKE_CURRENT_BINARY_DIR}/test-results/gemma_emmvgwas_compat_smoke"
    LABELS "compat;gemma;emmvgwas;smoke"
    TIMEOUT 120
)

add_custom_target(check_gemma_compat
    COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_CURRENT_BINARY_DIR} -R gemma_emmvgwas_compat_smoke --output-on-failure
    DEPENDS cosmicscanner
    USES_TERMINAL
)

add_test(
    NAME gemma_bslmm_reference_validation
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_gemma_bslmm_reference_validation.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
set_tests_properties(gemma_bslmm_reference_validation PROPERTIES
    ENVIRONMENT "COSMICSCANNER_ROOT=${CMAKE_CURRENT_SOURCE_DIR};COSMICSCANNER_BIN=$<TARGET_FILE:cosmicscanner>;COSMIC_BSLMM_REF_OUT=${CMAKE_CURRENT_BINARY_DIR}/test-results/gemma_bslmm_reference_validation"
    LABELS "reference;gemma;bslmm"
    SKIP_RETURN_CODE 2
    TIMEOUT 180
)

add_custom_target(check_gemma_reference
    COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_CURRENT_BINARY_DIR} -R gemma_bslmm_reference_validation --output-on-failure
    DEPENDS cosmicscanner
    USES_TERMINAL
)

add_test(
    NAME gemma_multivariate_reference_validation
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_gemma_multivariate_validation.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
set_tests_properties(gemma_multivariate_reference_validation PROPERTIES
    ENVIRONMENT "COSMICSCANNER_ROOT=${CMAKE_CURRENT_SOURCE_DIR};COSMICSCANNER_BIN=$<TARGET_FILE:cosmicscanner>"
    LABELS "reference;gemma;emmvgwas;multivariate"
    SKIP_RETURN_CODE 2
    TIMEOUT 300
)

add_custom_target(check_gemma_mv_reference
    COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_CURRENT_BINARY_DIR} -R gemma_multivariate_reference_validation --output-on-failure
    DEPENDS cosmicscanner
    USES_TERMINAL
)

add_custom_target(check_gemma_all
    COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_CURRENT_BINARY_DIR} -L gemma --output-on-failure
    DEPENDS cosmicscanner
    USES_TERMINAL
)

add_test(
    NAME eb_gwas_suite_quick
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_eb_gwas_suite.sh --quick
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
set_tests_properties(eb_gwas_suite_quick PROPERTIES
    ENVIRONMENT "COSMICSCANNER_BIN=$<TARGET_FILE:cosmicscanner>"
    LABELS "eb;gwas;smoke;validation"
    TIMEOUT 180
)

add_custom_target(check_eb_gwas
    COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_CURRENT_BINARY_DIR} -R eb_gwas_suite_quick --output-on-failure
    DEPENDS cosmicscanner
    USES_TERMINAL
)

add_custom_target(check_eb_gwas_full
    COMMAND ${CMAKE_COMMAND} -E env COSMICSCANNER_BIN=$<TARGET_FILE:cosmicscanner> bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_eb_gwas_suite.sh --full
    DEPENDS cosmicscanner
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    USES_TERMINAL
)

add_test(
    NAME ssgwas_a22_smoke
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_test_singlestep.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
set_tests_properties(ssgwas_a22_smoke PROPERTIES
    ENVIRONMENT "COSMICSCANNER_BIN=$<TARGET_FILE:cosmicscanner>"
    LABELS "ssgwas;gwas;smoke;validation"
    TIMEOUT 120
)

add_custom_target(check_ssgwas
    COMMAND ${CMAKE_CTEST_COMMAND} --test-dir ${CMAKE_CURRENT_BINARY_DIR} -R ssgwas_a22_smoke --output-on-failure
    DEPENDS cosmicscanner
    USES_TERMINAL
)

add_executable(test_sparse_reml tests/test_sparse_reml.cpp)
target_link_libraries(test_sparse_reml PRIVATE CosmicCommon)
add_test(NAME sparse_reml_unit COMMAND test_sparse_reml)
set_tests_properties(sparse_reml_unit PROPERTIES LABELS "unit;vce" TIMEOUT 60)

add_executable(test_multi_vc tests/test_multi_vc.cpp)
target_link_libraries(test_multi_vc PRIVATE CosmicCommon)
add_test(NAME multi_vc_unit COMMAND test_multi_vc)
set_tests_properties(multi_vc_unit PROPERTIES LABELS "unit;vce" TIMEOUT 60)

add_executable(test_spa_accuracy tests/test_spa_accuracy.cpp src/spa_test.cpp)
target_link_libraries(test_spa_accuracy CosmicCommon)
target_include_directories(test_spa_accuracy PRIVATE src)
add_test(NAME spa_accuracy_unit COMMAND test_spa_accuracy)
set_tests_properties(spa_accuracy_unit PROPERTIES LABELS "unit;spa" TIMEOUT 60)

add_executable(test_regenie tests/test_regenie.cpp src/ridge_regenie.cpp)
target_link_libraries(test_regenie CosmicCommon)
target_include_directories(test_regenie PRIVATE src)
add_test(NAME regenie_unit COMMAND test_regenie)
set_tests_properties(regenie_unit PROPERTIES LABELS "unit;regenie" TIMEOUT 60)

add_executable(test_input_adaptation_units
    tests/test_input_adaptation_units.cpp
    src/input/input_common.cpp
    src/input/input_source_reader.cpp
    src/input/time_normalizer.cpp
    src/input/wide_to_long.cpp
)
target_include_directories(test_input_adaptation_units PRIVATE src)
add_test(NAME input_adaptation_unit COMMAND test_input_adaptation_units)
set_tests_properties(input_adaptation_unit PROPERTIES LABELS "unit;input" TIMEOUT 60)

add_executable(eval_repeatability_reml tools/eval_repeatability_reml.cpp)
target_link_libraries(eval_repeatability_reml PRIVATE CosmicCommon)
target_include_directories(eval_repeatability_reml PRIVATE
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)
endif()

install(TARGETS cosmicscanner RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
