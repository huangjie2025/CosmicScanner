cmake_minimum_required(VERSION 3.15)
project(CosmicScanner VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(WIN32)
    add_compile_definitions(NOMINMAX WIN32_LEAN_AND_MEAN)
endif()
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
file(GLOB ASSOCIATION_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/association/*.cpp")
file(GLOB ASSOCIATION_MULTILOCUS_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/association/multilocus/*.cpp")
file(GLOB ASSOCIATION_GENERALIZED_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/association/generalized/*.cpp")
file(GLOB ASSOCIATION_COJO_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/association/cojo/*.cpp")
file(GLOB ASSOCIATION_META_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/association/meta/*.cpp")
file(GLOB ASSOCIATION_SUMSTATS_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/association/sumstats/*.cpp")
file(GLOB ASSOCIATION_FINEMAP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/association/fine_mapping/*.cpp")
file(GLOB VCF_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/libCosmicCore/src/vcf/*.cpp")
file(GLOB GENETICS_LD_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/genetics/ld/*.cpp")
file(GLOB GENETICS_QC_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/genetics/qc/*.cpp")
file(GLOB GENETICS_POP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/genetics/population/*.cpp")
list(APPEND COMMON_SOURCES ${MODEL_SOURCES})
list(APPEND COMMON_SOURCES ${INPUT_SOURCES})
list(APPEND COMMON_SOURCES ${ASSOCIATION_SOURCES})
list(APPEND COMMON_SOURCES ${ASSOCIATION_MULTILOCUS_SOURCES})
list(APPEND COMMON_SOURCES ${ASSOCIATION_GENERALIZED_SOURCES})
list(APPEND COMMON_SOURCES ${ASSOCIATION_COJO_SOURCES})
list(APPEND COMMON_SOURCES ${ASSOCIATION_META_SOURCES})
list(APPEND COMMON_SOURCES ${ASSOCIATION_SUMSTATS_SOURCES})
list(APPEND COMMON_SOURCES ${ASSOCIATION_FINEMAP_SOURCES})
list(APPEND COMMON_SOURCES ${GENETICS_LD_SOURCES})
list(APPEND COMMON_SOURCES ${GENETICS_QC_SOURCES})
list(APPEND COMMON_SOURCES ${GENETICS_POP_SOURCES})
list(APPEND COMMON_SOURCES ${VCF_SOURCES})
list(APPEND COMMON_SOURCES src/models/rare_variant/rare_variant_extensions.cpp)
# 过滤掉带有 main 函数的测试或入口文件，防止链接冲突
list(REMOVE_ITEM COMMON_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
list(REMOVE_ITEM COMMON_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/gwas_app.cpp")

add_executable(cosmicscanner src/main.cpp src/gwas_app.cpp ${COMMON_SOURCES})
if(WIN32)
    enable_language(RC)
    target_sources(cosmicscanner PRIVATE resources/windows_manifest.rc)
    target_include_directories(cosmicscanner PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/resources)
endif()

target_include_directories(cosmicscanner PRIVATE 
    ${CMAKE_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/src/models
    ${CMAKE_SOURCE_DIR}/src/association
    ${CMAKE_SOURCE_DIR}/src/genetics
    ${CMAKE_SOURCE_DIR}/src/genetics/ld
    ${CMAKE_SOURCE_DIR}/src/genetics/qc
    ${CMAKE_SOURCE_DIR}/src/genetics/population
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)

target_link_libraries(cosmicscanner PRIVATE CosmicCommon OpenMP::OpenMP_CXX)
target_compile_definitions(cosmicscanner PRIVATE COSMICSCANNER_VERSION="${PROJECT_VERSION}")

if(COSMICSCANNER_BUILD_TESTS)
enable_testing()
find_package(Python3 COMPONENTS Interpreter REQUIRED)
add_test(NAME cli_regressions COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_cli_regressions.py
    --bin $<TARGET_FILE:cosmicscanner>)
set_tests_properties(cli_regressions PROPERTIES LABELS "cli;regression;integration" TIMEOUT 180
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
add_test(NAME review_regressions COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_review_regressions.py
    --bin $<TARGET_FILE:cosmicscanner>)
set_tests_properties(review_regressions PROPERTIES LABELS "review;regression;integration" TIMEOUT 180
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
add_executable(profile_kinship_kernel tests/profile_kinship_kernel.cpp src/kinship_operator.cpp)
target_include_directories(profile_kinship_kernel PRIVATE src)
target_link_libraries(profile_kinship_kernel PRIVATE CosmicCommon OpenMP::OpenMP_CXX)
add_executable(test_kinship_batch tests/test_kinship_batch.cpp src/kinship_operator.cpp)
target_include_directories(test_kinship_batch PRIVATE src)
target_link_libraries(test_kinship_batch PRIVATE CosmicCommon OpenMP::OpenMP_CXX)
add_test(NAME kinship_batch_equivalence COMMAND test_kinship_batch ${CMAKE_CURRENT_BINARY_DIR}/test-results/kinship-batch)
set_tests_properties(kinship_batch_equivalence PROPERTIES TIMEOUT 120
    ENVIRONMENT "OMP_NUM_THREADS=4;MKL_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1")
add_test(NAME r27_release_closure COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_r27_release_closure.py
    --bin $<TARGET_FILE:cosmicscanner> --out ${CMAKE_CURRENT_BINARY_DIR}/test-results/r27)
set_tests_properties(r27_release_closure PROPERTIES LABELS "release;integration" TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
add_executable(test_population_ld tests/test_population_ld.cpp ${GENETICS_POP_SOURCES} ${GENETICS_LD_SOURCES})
target_include_directories(test_population_ld PRIVATE src)
target_link_libraries(test_population_ld PRIVATE CosmicCommon)
add_test(NAME population_ld_unit COMMAND test_population_ld)
add_executable(test_qc_pipeline tests/test_qc_pipeline.cpp ${GENETICS_QC_SOURCES})
target_include_directories(test_qc_pipeline PRIVATE src)
target_link_libraries(test_qc_pipeline PRIVATE CosmicCommon)
add_test(NAME qc_pipeline_unit COMMAND test_qc_pipeline ${CMAKE_CURRENT_BINARY_DIR}/test-results/qc-unit)
add_executable(test_association_interface tests/test_association_interface.cpp ${ASSOCIATION_SOURCES} ${ASSOCIATION_MULTILOCUS_SOURCES} ${ASSOCIATION_GENERALIZED_SOURCES} ${GENETICS_LD_SOURCES} src/spa_test.cpp)
target_include_directories(test_association_interface PRIVATE src)
target_link_libraries(test_association_interface PRIVATE CosmicCommon)
add_test(NAME association_interface_unit COMMAND test_association_interface)
add_executable(test_sumstats_harmonization tests/test_sumstats_harmonization.cpp ${ASSOCIATION_SUMSTATS_SOURCES})
target_include_directories(test_sumstats_harmonization PRIVATE src)
target_link_libraries(test_sumstats_harmonization PRIVATE CosmicCommon)
add_test(NAME sumstats_harmonization_unit COMMAND test_sumstats_harmonization)
set_tests_properties(sumstats_harmonization_unit PROPERTIES LABELS "sumstats;harmonization;unit" TIMEOUT 30)
add_executable(test_rare_variant_extensions tests/test_rare_variant_extensions.cpp src/rare_variant_test.cpp src/models/rare_variant/rare_variant_extensions.cpp)
target_include_directories(test_rare_variant_extensions PRIVATE src)
target_link_libraries(test_rare_variant_extensions PRIVATE CosmicCommon)
add_test(NAME rare_variant_extensions_unit COMMAND test_rare_variant_extensions)
add_test(NAME r8_b2_binary_rare_integration COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_r8_b2_rare.py
    --bin $<TARGET_FILE:cosmicscanner> --out ${CMAKE_CURRENT_BINARY_DIR}/test-results/r8-b2)
set_tests_properties(r8_b2_binary_rare_integration PROPERTIES LABELS "association;r8-b2;integration" TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
add_test(NAME r8_b2_rare_vcf COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_r8_b2_rare.py
    --bin $<TARGET_FILE:cosmicscanner> --out ${CMAKE_CURRENT_BINARY_DIR}/test-results/r8-b2-rare)
set_tests_properties(r8_b2_rare_vcf PROPERTIES LABELS "association;r8-b2;integration" TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
set_tests_properties(population_ld_unit qc_pipeline_unit association_interface_unit rare_variant_extensions_unit
    PROPERTIES LABELS "association-r2;unit" TIMEOUT 120)
foreach(case pca qc_filter ld_prune ld_clump skat_o acat_v genetic_model)
    add_test(NAME ${case}_smoke COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_association_r2.py
        --bin $<TARGET_FILE:cosmicscanner> --case ${case}
        --out ${CMAKE_CURRENT_BINARY_DIR}/test-results/association-r2/${case})
    set_tests_properties(${case}_smoke PROPERTIES LABELS "association-r2;smoke" TIMEOUT 180
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
endforeach()
add_test(NAME ad_gwas_p2b_validation COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ad_gwas_validation.py
    --scanner $<TARGET_FILE:cosmicscanner>
    --workdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/ad-gwas-p2b)
set_tests_properties(ad_gwas_p2b_validation PROPERTIES LABELS "association;ad;integration" TIMEOUT 600
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
add_test(NAME ad_gwas_cli_guard COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ad_cli_guard_tests.py
    --scanner $<TARGET_FILE:cosmicscanner>
    --workdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/ad-gwas-cli-guard)
set_tests_properties(ad_gwas_cli_guard PROPERTIES LABELS "association;ad;cli" TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
add_executable(test_multilocus_engine tests/test_multilocus_engine.cpp ${ASSOCIATION_SOURCES} ${ASSOCIATION_MULTILOCUS_SOURCES} ${ASSOCIATION_GENERALIZED_SOURCES} ${GENETICS_LD_SOURCES} src/spa_test.cpp)
target_include_directories(test_multilocus_engine PRIVATE src)
target_link_libraries(test_multilocus_engine PRIVATE CosmicCommon)
add_test(NAME multilocus_engine_unit COMMAND test_multilocus_engine)
set_tests_properties(multilocus_engine_unit PROPERTIES LABELS "association-r3;association-r4a;unit" TIMEOUT 120)
add_executable(test_firth_logistic tests/test_firth_logistic.cpp ${ASSOCIATION_SOURCES} ${ASSOCIATION_MULTILOCUS_SOURCES} ${ASSOCIATION_GENERALIZED_SOURCES} ${GENETICS_LD_SOURCES} src/spa_test.cpp)
target_include_directories(test_firth_logistic PRIVATE src)
target_link_libraries(test_firth_logistic PRIVATE CosmicCommon)
add_test(NAME firth_logistic_unit COMMAND test_firth_logistic)
set_tests_properties(firth_logistic_unit PROPERTIES LABELS "association-r4a;unit" TIMEOUT 120)
add_executable(test_saige_style tests/test_saige_style.cpp ${ASSOCIATION_SOURCES} ${ASSOCIATION_MULTILOCUS_SOURCES} ${ASSOCIATION_GENERALIZED_SOURCES} ${GENETICS_LD_SOURCES} src/spa_test.cpp)
target_include_directories(test_saige_style PRIVATE src)
target_link_libraries(test_saige_style PRIVATE CosmicCommon)
add_test(NAME saige_style_unit COMMAND test_saige_style)
set_tests_properties(saige_style_unit PROPERTIES LABELS "association-r5a;unit" TIMEOUT 120)
add_executable(test_logistic tests/test_logistic.cpp ${ASSOCIATION_SOURCES} ${ASSOCIATION_MULTILOCUS_SOURCES} ${ASSOCIATION_GENERALIZED_SOURCES} ${GENETICS_LD_SOURCES} src/spa_test.cpp)
target_include_directories(test_logistic PRIVATE src)
target_link_libraries(test_logistic PRIVATE CosmicCommon)
add_test(NAME logistic_unit COMMAND test_logistic)
set_tests_properties(logistic_unit PROPERTIES LABELS "association;r8-b2;unit" TIMEOUT 120)
add_executable(test_genotype_encoding tests/test_genotype_encoding.cpp)
target_include_directories(test_genotype_encoding PRIVATE src)
target_link_libraries(test_genotype_encoding PRIVATE CosmicCommon)
add_test(NAME genotype_encoding_unit COMMAND test_genotype_encoding)
set_tests_properties(genotype_encoding_unit PROPERTIES LABELS "association;r8-b3;unit" TIMEOUT 30)
add_test(NAME low_rank_lmm_validation COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_low_rank_lmm_validation.py
    --bin $<TARGET_FILE:cosmicscanner> --out ${CMAKE_CURRENT_BINARY_DIR}/test-results/r8-b3-low-rank)
set_tests_properties(low_rank_lmm_validation PROPERTIES LABELS "association;r8-b3;integration" TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
add_executable(test_new_readers tests/test_new_readers.cpp)
target_link_libraries(test_new_readers PRIVATE CosmicCommon)
add_test(NAME new_readers_unit COMMAND test_new_readers ${CMAKE_CURRENT_BINARY_DIR}/test-results/new-readers)
set_tests_properties(new_readers_unit PROPERTIES LABELS "reader;r6-c;unit" TIMEOUT 30)
add_test(NAME new_format_validation COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_new_format_validation.py
    --bin $<TARGET_FILE:cosmicscanner> --out ${CMAKE_CURRENT_BINARY_DIR}/test-results/r6-c-new-formats)
set_tests_properties(new_format_validation PROPERTIES LABELS "reader;r6-c;integration" TIMEOUT 180
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
add_executable(test_cojo tests/test_cojo.cpp ${ASSOCIATION_SOURCES} ${ASSOCIATION_MULTILOCUS_SOURCES} ${ASSOCIATION_GENERALIZED_SOURCES} ${ASSOCIATION_COJO_SOURCES} ${ASSOCIATION_META_SOURCES} ${ASSOCIATION_SUMSTATS_SOURCES} ${GENETICS_LD_SOURCES} ${VCF_SOURCES} src/spa_test.cpp)
target_include_directories(test_cojo PRIVATE src)
target_link_libraries(test_cojo PRIVATE CosmicCommon)
add_test(NAME cojo_unit COMMAND test_cojo)
set_tests_properties(cojo_unit PROPERTIES LABELS "association-r5b;unit" TIMEOUT 120)
add_test(NAME cojo_cli_validation COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_cojo_validation.py
    --bin $<TARGET_FILE:cosmicscanner> --out ${CMAKE_CURRENT_BINARY_DIR}/test-results/r8-b1-cojo)
set_tests_properties(cojo_cli_validation PROPERTIES LABELS "association-r5b;cojo;r8-b1;integration" TIMEOUT 120
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
add_executable(test_fine_mapping tests/test_fine_mapping.cpp ${ASSOCIATION_FINEMAP_SOURCES} ${ASSOCIATION_SUMSTATS_SOURCES})
target_include_directories(test_fine_mapping PRIVATE src)
target_link_libraries(test_fine_mapping PRIVATE CosmicCommon)
add_test(NAME fine_mapping_unit COMMAND test_fine_mapping)
set_tests_properties(fine_mapping_unit PROPERTIES LABELS "fine-mapping;r7-a;unit" TIMEOUT 30)
add_test(NAME fine_mapping_cli_validation COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_fine_mapping_validation.py
    --bin $<TARGET_FILE:cosmicscanner> --out ${CMAKE_CURRENT_BINARY_DIR}/test-results/r7-a-fine-mapping)
set_tests_properties(fine_mapping_cli_validation PROPERTIES LABELS "fine-mapping;r7-a;integration" TIMEOUT 120
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
add_executable(test_meta_analysis tests/test_meta_analysis.cpp ${ASSOCIATION_SOURCES} ${ASSOCIATION_MULTILOCUS_SOURCES} ${ASSOCIATION_GENERALIZED_SOURCES} ${ASSOCIATION_COJO_SOURCES} ${ASSOCIATION_META_SOURCES} ${ASSOCIATION_SUMSTATS_SOURCES} ${GENETICS_LD_SOURCES} ${VCF_SOURCES} src/spa_test.cpp)
target_include_directories(test_meta_analysis PRIVATE src)
target_link_libraries(test_meta_analysis PRIVATE CosmicCommon)
add_test(NAME meta_analysis_unit COMMAND test_meta_analysis)
set_tests_properties(meta_analysis_unit PROPERTIES LABELS "association-r5b;unit" TIMEOUT 120)
foreach(case farmcpu blink conditional)
    add_test(NAME ${case}_smoke_test COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_association_r3.py
        --bin $<TARGET_FILE:cosmicscanner> --case ${case}
        --out ${CMAKE_CURRENT_BINARY_DIR}/test-results/association-r3/${case})
    set_tests_properties(${case}_smoke_test PROPERTIES LABELS "association-r3;smoke" TIMEOUT 180
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
endforeach()
add_custom_target(benchmark_association_r3 COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_association_r3.py
    --bin $<TARGET_FILE:cosmicscanner> --case benchmark
    --out ${CMAKE_CURRENT_BINARY_DIR}/test-results/association-r3/benchmark
    DEPENDS cosmicscanner USES_TERMINAL)
foreach(case firth conditional_joint)
    add_test(NAME ${case}_smoke_test COMMAND ${Python3_EXECUTABLE}
        ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_association_r4a.py
        --bin $<TARGET_FILE:cosmicscanner> --case ${case}
        --out ${CMAKE_CURRENT_BINARY_DIR}/test-results/association-r4a/${case})
    set_tests_properties(${case}_smoke_test PROPERTIES LABELS "association-r4a;smoke" TIMEOUT 180
        ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
endforeach()
add_custom_target(benchmark_association_r4a COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_association_r4a.py
    --bin $<TARGET_FILE:cosmicscanner> --case benchmark
    --out ${CMAKE_CURRENT_BINARY_DIR}/test-results/association-r4a/benchmark
    DEPENDS cosmicscanner USES_TERMINAL)
add_test(NAME saige_style_smoke_test COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_association_r5a.py
    --bin $<TARGET_FILE:cosmicscanner> --case saige_style
    --out ${CMAKE_CURRENT_BINARY_DIR}/test-results/association-r5a/saige_style)
set_tests_properties(saige_style_smoke_test PROPERTIES LABELS "association-r5a;smoke" TIMEOUT 180
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
add_custom_target(benchmark_association_r5a COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_association_r5a.py
    --bin $<TARGET_FILE:cosmicscanner> --case benchmark
    --out ${CMAKE_CURRENT_BINARY_DIR}/test-results/association-r5a/benchmark
    DEPENDS cosmicscanner USES_TERMINAL)
add_custom_target(benchmark_association_r2 COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_association_r2.py
    --bin $<TARGET_FILE:cosmicscanner> --case benchmark
    --out ${CMAKE_CURRENT_BINARY_DIR}/test-results/association-r2/benchmark
    DEPENDS cosmicscanner USES_TERMINAL)

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
    NAME ssgwas_backsolve_validation
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ssgwas_backsolve.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
set_tests_properties(ssgwas_backsolve_validation PROPERTIES
    ENVIRONMENT "COSMICSCANNER_ROOT=${CMAKE_CURRENT_SOURCE_DIR};COSMICSCANNER_BIN=$<TARGET_FILE:cosmicscanner>;COSMIC_BACKSOLVE_OUT=${CMAKE_CURRENT_BINARY_DIR}/test-results/ssgwas_backsolve"
    LABELS "ssgwas;backsolve;analytic;simulation"
    TIMEOUT 900
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

add_test(
    NAME analysis_dataset_lmm_migration
    COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_analysis_dataset_lmm_migration.sh
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)
set_tests_properties(analysis_dataset_lmm_migration PROPERTIES
    ENVIRONMENT "COSMICSCANNER_ROOT=${CMAKE_CURRENT_SOURCE_DIR};COSMICSCANNER_BIN=$<TARGET_FILE:cosmicscanner>;COSMIC_ANALYSIS_DATASET_LMM_OUT=${CMAKE_CURRENT_BINARY_DIR}/test-results/analysis_dataset_lmm"
    LABELS "analysis_dataset;gwas;lmm;unit"
    TIMEOUT 120
)

add_executable(test_sparse_reml tests/test_sparse_reml.cpp)
target_link_libraries(test_sparse_reml PRIVATE CosmicCommon)
add_test(NAME sparse_reml_unit COMMAND test_sparse_reml)
set_tests_properties(sparse_reml_unit PROPERTIES LABELS "unit;vce" TIMEOUT 60)

add_executable(test_multi_vc tests/test_multi_vc.cpp)
target_link_libraries(test_multi_vc PRIVATE CosmicCommon)
add_test(NAME multi_vc_unit COMMAND test_multi_vc)
set_tests_properties(multi_vc_unit PROPERTIES LABELS "unit;vce" TIMEOUT 60)

add_executable(test_mvlmm_numerics tests/test_mvlmm_numerics.cpp)
target_link_libraries(test_mvlmm_numerics PRIVATE CosmicCommon)
target_include_directories(test_mvlmm_numerics PRIVATE src)
add_test(NAME mvlmm_numerics_unit COMMAND test_mvlmm_numerics)
set_tests_properties(mvlmm_numerics_unit PROPERTIES LABELS "unit;mvlmm;numerics" TIMEOUT 60)

add_executable(test_rrm_basis_identity tests/test_rrm_basis_identity.cpp)
target_link_libraries(test_rrm_basis_identity PRIVATE CosmicCommon)
add_test(NAME rrm_basis_identity_unit COMMAND test_rrm_basis_identity)
set_tests_properties(rrm_basis_identity_unit PROPERTIES LABELS "unit;rrm" TIMEOUT 60)

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

add_executable(test_apy_ginv
    tests/test_apy_ginv.cpp
    src/models/apy_ginv.cpp
    src/models/gcc_solver.cpp
    src/models/gcc_solver_pcg.cpp
    src/models/gcc_preconditioner.cpp
)
target_link_libraries(test_apy_ginv PRIVATE CosmicCommon)
target_include_directories(test_apy_ginv PRIVATE
    src
    src/models
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)
add_test(NAME apy_ginv_unit COMMAND test_apy_ginv)
set_tests_properties(apy_ginv_unit PROPERTIES LABELS "unit;ssgwas;apy" TIMEOUT 120)

add_executable(test_e1_operator tests/test_e1_operator.cpp ${VCF_SOURCES})
target_link_libraries(test_e1_operator PRIVATE CosmicCommon)
target_include_directories(test_e1_operator PRIVATE
    src
    src/models
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)
add_test(NAME e1_operator_unit COMMAND test_e1_operator)
set_tests_properties(e1_operator_unit PROPERTIES LABELS "unit;ssgwas;apy;e1" TIMEOUT 300)

add_executable(test_e2_stream_parity
    tests/test_e2_stream_parity.cpp
    src/models/block_streamed_grm.cpp
    src/models/apy_ginv.cpp
    src/models/gcc_solver.cpp
    src/models/gcc_solver_pcg.cpp
    src/models/gcc_preconditioner.cpp
    ${VCF_SOURCES}
)
target_link_libraries(test_e2_stream_parity PRIVATE CosmicCommon)
target_include_directories(test_e2_stream_parity PRIVATE
    src
    src/models
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)
add_test(NAME e2_stream_parity_unit COMMAND test_e2_stream_parity)
set_tests_properties(e2_stream_parity_unit PROPERTIES LABELS "unit;ssgwas;apy;e2" TIMEOUT 300)

add_executable(test_gcc_solver_parity
    tests/test_gcc_solver_parity.cpp
    src/models/apy_ginv.cpp
    src/models/gcc_solver.cpp
    src/models/gcc_solver_pcg.cpp
    src/models/gcc_preconditioner.cpp
    ${VCF_SOURCES}
)
target_link_libraries(test_gcc_solver_parity PRIVATE CosmicCommon)
target_include_directories(test_gcc_solver_parity PRIVATE
    src
    src/models
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)
add_test(NAME gcc_solver_parity_unit COMMAND test_gcc_solver_parity)
set_tests_properties(gcc_solver_parity_unit PROPERTIES LABELS "unit;ssgwas;apy;e3a" TIMEOUT 300)

add_executable(test_gcc_pcg_solver
    tests/test_gcc_pcg_solver.cpp
    src/models/gcc_solver.cpp
    src/models/gcc_solver_pcg.cpp
    src/models/gcc_preconditioner.cpp
    ${VCF_SOURCES}
)
target_link_libraries(test_gcc_pcg_solver PRIVATE CosmicCommon)
target_include_directories(test_gcc_pcg_solver PRIVATE
    src
    src/models
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)
add_test(NAME gcc_pcg_solver_unit COMMAND test_gcc_pcg_solver)
set_tests_properties(gcc_pcg_solver_unit PROPERTIES LABELS "unit;ssgwas;apy;e3b1" TIMEOUT 300)

add_executable(test_apy_operator
    tests/test_apy_operator.cpp
    src/models/apy_ginv.cpp
    src/models/apy_inverse_operator.cpp
    src/models/cross_genomic_operator.cpp
    src/models/gcc_solver.cpp
    src/models/gcc_solver_pcg.cpp
    src/models/gcc_preconditioner.cpp
    ${VCF_SOURCES}
)
target_link_libraries(test_apy_operator PRIVATE CosmicCommon OpenMP::OpenMP_CXX)
target_include_directories(test_apy_operator PRIVATE
    src
    src/models
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)
add_test(NAME apy_operator_unit COMMAND test_apy_operator)
set_tests_properties(apy_operator_unit PROPERTIES LABELS "unit;ssgwas;apy;e3c2" TIMEOUT 300)
add_test(NAME apy_operator_downstream_parity COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ssgwas_apy_operator_parity.py
    --bin $<TARGET_FILE:cosmicscanner>
    --outdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/apy-operator-parity)
set_tests_properties(apy_operator_downstream_parity PROPERTIES
    LABELS "integration;ssgwas;apy;e3c2" TIMEOUT 300
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

add_executable(test_gcn_genotype_streaming
    tests/test_gcn_genotype_streaming.cpp
    src/models/apy_ginv.cpp
    src/models/apy_inverse_operator.cpp
    src/models/cross_genomic_operator.cpp
    src/models/block_streamed_grm.cpp
    src/models/gcc_solver.cpp
    src/models/gcc_solver_pcg.cpp
    src/models/gcc_preconditioner.cpp
    ${VCF_SOURCES}
)
target_link_libraries(test_gcn_genotype_streaming PRIVATE CosmicCommon OpenMP::OpenMP_CXX)
target_include_directories(test_gcn_genotype_streaming PRIVATE
    src
    src/models
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)
add_test(NAME gcn_genotype_streaming_unit COMMAND test_gcn_genotype_streaming)
set_tests_properties(gcn_genotype_streaming_unit PROPERTIES
    LABELS "unit;ssgwas;apy;e3c3" TIMEOUT 300)
add_test(NAME gcn_genotype_streaming_frozen COMMAND test_gcn_genotype_streaming
    ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/ssgwas_blupf90_parity/dataset/sim)
set_tests_properties(gcn_genotype_streaming_frozen PROPERTIES
    LABELS "integration;ssgwas;apy;e3c3" TIMEOUT 300)

add_executable(test_apy_mme_operator
    tests/test_apy_mme_operator.cpp
    src/models/apy_ginv.cpp
    src/models/apy_inverse_operator.cpp
    src/models/cross_genomic_operator.cpp
    src/models/gcc_solver.cpp
    src/models/gcc_solver_pcg.cpp
    src/models/gcc_preconditioner.cpp
    ${VCF_SOURCES}
)
target_link_libraries(test_apy_mme_operator PRIVATE CosmicCommon OpenMP::OpenMP_CXX)
target_include_directories(test_apy_mme_operator PRIVATE
    src
    src/models
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)
add_test(NAME apy_mme_operator_unit COMMAND test_apy_mme_operator)
set_tests_properties(apy_mme_operator_unit PROPERTIES
    LABELS "unit;ssgwas;apy;e3d1" TIMEOUT 300)
add_test(NAME apy_mme_operator_downstream_parity COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ssgwas_apy_mme_operator_parity.py
    --bin $<TARGET_FILE:cosmicscanner>
    --outdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/apy-mme-operator-parity)
set_tests_properties(apy_mme_operator_downstream_parity PROPERTIES
    LABELS "integration;ssgwas;apy;e3d1" TIMEOUT 600
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

add_test(NAME pev_scan_scalar_batched_parity COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ssgwas_pev_scan_parity.py
    --bin $<TARGET_FILE:cosmicscanner>
    --outdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/pev-scan-parity)
set_tests_properties(pev_scan_scalar_batched_parity PROPERTIES
    LABELS "integration;ssgwas;apy;e5" TIMEOUT 1200
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

add_executable(test_a22_hinv_operator
    tests/test_a22_hinv_operator.cpp
    src/models/a22_inverse_operator.cpp
    src/models/h_inverse_operator.cpp
    src/models/pedigree_relationship.cpp
)
target_link_libraries(test_a22_hinv_operator PRIVATE CosmicCommon OpenMP::OpenMP_CXX)
target_include_directories(test_a22_hinv_operator PRIVATE
    src
    src/models
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)
add_test(NAME a22_hinv_operator_unit COMMAND test_a22_hinv_operator)
set_tests_properties(a22_hinv_operator_unit PROPERTIES
    LABELS "unit;ssgwas;a22;hinv;e4c" TIMEOUT 300)
add_test(NAME a22_hinv_operator_downstream_parity COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ssgwas_hinv_operator_parity.py
    --bin $<TARGET_FILE:cosmicscanner>
    --outdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/a22-hinv-operator-parity)
set_tests_properties(a22_hinv_operator_downstream_parity PROPERTIES
    LABELS "integration;ssgwas;a22;hinv;e4c" TIMEOUT 600
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

add_executable(test_pev_operator
    tests/test_pev_operator.cpp
    src/models/apy_core_pev_operator.cpp
    src/models/pev_operator.cpp
    src/models/apy_ginv.cpp
    src/models/apy_inverse_operator.cpp
    src/models/cross_genomic_operator.cpp
    src/models/gcc_solver.cpp
    src/models/gcc_solver_pcg.cpp
    src/models/gcc_preconditioner.cpp
    ${VCF_SOURCES}
)
target_link_libraries(test_pev_operator PRIVATE CosmicCommon OpenMP::OpenMP_CXX)
target_include_directories(test_pev_operator PRIVATE
    src
    src/models
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)
add_test(NAME pev_operator_unit COMMAND test_pev_operator)
set_tests_properties(pev_operator_unit PROPERTIES
    LABELS "unit;ssgwas;apy;e3d2b;e3d2c" TIMEOUT 300)
add_test(NAME pev_operator_frozen COMMAND test_pev_operator
    ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/ssgwas_blupf90_parity/dataset/sim)
set_tests_properties(pev_operator_frozen PROPERTIES
    LABELS "integration;ssgwas;apy;e3d2b;e3d2c" TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
add_test(NAME apy_pev_operator_downstream_parity COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ssgwas_apy_pev_operator_parity.py
    --bin $<TARGET_FILE:cosmicscanner>
    --outdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/apy-pev-operator-parity)
set_tests_properties(apy_pev_operator_downstream_parity PROPERTIES
    LABELS "integration;ssgwas;apy;e3d2c" TIMEOUT 900
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

add_executable(test_dominance_scan tests/test_dominance_scan.cpp src/dominance_engine.cpp ${VCF_SOURCES})
target_include_directories(test_dominance_scan PRIVATE src ${CMAKE_SOURCE_DIR}/libCosmicCore/src)
target_link_libraries(test_dominance_scan PRIVATE CosmicCommon)
add_test(NAME dominance_scan_unit COMMAND test_dominance_scan)
set_tests_properties(dominance_scan_unit PROPERTIES LABELS "unit;dominance;scan" TIMEOUT 120)
add_test(NAME dominance_gwas_validation COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_dominance_gwas_validation.py
    --scanner $<TARGET_FILE:cosmicscanner>
    --workdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/dominance-gwas)
set_tests_properties(dominance_gwas_validation PROPERTIES
    LABELS "integration;dominance;validation" TIMEOUT 900
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")
add_test(NAME pgen_dominance_ad_validation COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_pgen_dominance_ad_validation.py
    --scanner $<TARGET_FILE:cosmicscanner>
    --workdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/pgen-dominance-ad
    --n 500 --m 600)
set_tests_properties(pgen_dominance_ad_validation PROPERTIES
    LABELS "integration;dominance;ad;pgen;e2.1;validation" TIMEOUT 2400
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

add_executable(test_bgen_probability_backend
    tests/test_bgen_probability_backend.cpp)
target_include_directories(test_bgen_probability_backend PRIVATE
    src ${CMAKE_SOURCE_DIR}/libCosmicCore/src ${EIGEN3_INCLUDE_DIR})
target_link_libraries(test_bgen_probability_backend PRIVATE CosmicCommon)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test_bgen_probability_backend PRIVATE OpenMP::OpenMP_CXX)
endif()
add_test(NAME bgen_probability_backend_unit COMMAND test_bgen_probability_backend)
set_tests_properties(bgen_probability_backend_unit PROPERTIES
    LABELS "unit;reader;bgen;e2.4" TIMEOUT 120)

add_test(NAME bgen_probability_cross_format COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_bgen_probability_cross_format.py
    --scanner $<TARGET_FILE:cosmicscanner>
    --workdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/bgen-prob-cross-format
    --n 400 --m 200)
set_tests_properties(bgen_probability_cross_format PROPERTIES
    LABELS "integration;reader;bgen;dominance;e2.6;validation" TIMEOUT 1200
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

add_test(NAME bgen_dominance_uncertainty_validation COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_bgen_dominance_uncertainty_validation.py
    --scanner $<TARGET_FILE:cosmicscanner>
    --workdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/bgen-dominance-uncertainty)
set_tests_properties(bgen_dominance_uncertainty_validation PROPERTIES
    LABELS "integration;bgen;dominance;e2.6;calibration" TIMEOUT 3600
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

add_test(NAME bgen_dominance_production_gate_smoke COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_bgen_dominance_production_gate.py
    --scanner $<TARGET_FILE:cosmicscanner>
    --workdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/bgen-dominance-production-smoke
    --n 200 --m 240 --reps 4 --workers 2 --smoke)
set_tests_properties(bgen_dominance_production_gate_smoke PROPERTIES
    LABELS "integration;bgen;dominance;e2.7;production" TIMEOUT 1200
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

add_test(NAME bgen_dominance_qc_validation COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_bgen_dominance_qc_validation.py
    --scanner $<TARGET_FILE:cosmicscanner>
    --workdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/bgen-dominance-qc)
set_tests_properties(bgen_dominance_qc_validation PROPERTIES
    LABELS "integration;bgen;dominance;e2.7;qc" TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

add_executable(test_bgen_ad_null tests/test_bgen_ad_null.cpp ${VCF_SOURCES})
target_link_libraries(test_bgen_ad_null PRIVATE CosmicCommon)
add_test(NAME bgen_ad_null_validation COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_bgen_ad_null_validation.py
    --driver $<TARGET_FILE:test_bgen_ad_null>
    --workdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/bgen-ad-null --legacy-control)
set_tests_properties(bgen_ad_null_validation PROPERTIES
    LABELS "integration;bgen;ad;null;e3.1" TIMEOUT 3600
    RESOURCE_LOCK "bgen_ad_null_fixture"
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

add_test(NAME bgen_ad_reml_robust_validation COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_bgen_ad_reml_robust_validation.py
    --driver $<TARGET_FILE:test_bgen_ad_null>
    --source ${CMAKE_CURRENT_BINARY_DIR}/test-results/bgen-ad-null
    --workdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/bgen-ad-reml-robust)
set_tests_properties(bgen_ad_reml_robust_validation PROPERTIES
    DEPENDS "bgen_ad_null_validation"
    LABELS "integration;bgen;ad;null;e3.1c;reml" TIMEOUT 3600
    RESOURCE_LOCK "bgen_ad_null_fixture"
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

add_executable(test_bgen_ad_scan tests/test_bgen_ad_scan.cpp src/bgen_ad_scan_engine.cpp src/ad_engine.cpp ${VCF_SOURCES})
target_include_directories(test_bgen_ad_scan PRIVATE src)
target_link_libraries(test_bgen_ad_scan PRIVATE CosmicCommon)
add_test(NAME bgen_ad_scan_engine_validation COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_bgen_ad_scan_validation.py
    --driver $<TARGET_FILE:test_bgen_ad_scan>
    --workdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/bgen-ad-scan)
set_tests_properties(bgen_ad_scan_engine_validation PROPERTIES
    LABELS "integration;bgen;ad;scan;e3.2a" TIMEOUT 7200
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

add_executable(test_bgen_ad_calibration tests/test_bgen_ad_calibration.cpp src/bgen_ad_scan_engine.cpp src/ad_engine.cpp ${VCF_SOURCES})
target_include_directories(test_bgen_ad_calibration PRIVATE src)
target_link_libraries(test_bgen_ad_calibration PRIVATE CosmicCommon OpenMP::OpenMP_CXX)

add_executable(test_bgen_ad_integration tests/test_bgen_ad_integration.cpp
    src/bgen_ad_gwas_scanner.cpp src/bgen_ad_scan_engine.cpp
    src/gwas_scanner.cpp src/ad_gwas_scanner.cpp src/ad_engine.cpp
    src/gwas_output.cpp src/gzip_writer.cpp src/lmm_engine.cpp src/mvlmm_engine.cpp
    src/matrix_free_vce.cpp src/kinship_operator.cpp src/spa_test.cpp
    src/association/sumstats/sumstats.cpp ${VCF_SOURCES})
target_include_directories(test_bgen_ad_integration PRIVATE src)
target_link_libraries(test_bgen_ad_integration PRIVATE CosmicCommon OpenMP::OpenMP_CXX)
add_test(NAME bgen_ad_scanner_integration COMMAND ${Python3_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_bgen_ad_integration_validation.py
    --harness $<TARGET_FILE:test_bgen_ad_integration> --scanner $<TARGET_FILE:cosmicscanner>
    --null-driver $<TARGET_FILE:test_bgen_ad_null>
    --workdir ${CMAKE_CURRENT_BINARY_DIR}/test-results/bgen-ad-integration --smoke)
set_tests_properties(bgen_ad_scanner_integration PROPERTIES
    LABELS "integration;bgen;ad;e3.2b;unified" TIMEOUT 10800
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1")

add_executable(test_cross_operator
    tests/test_cross_genomic_operator.cpp
    src/models/cross_genomic_operator.cpp
    src/models/block_streamed_grm.cpp
    src/models/apy_ginv.cpp
    src/models/apy_inverse_operator.cpp
    src/models/apy_core_pev_operator.cpp
    src/models/pev_operator.cpp
    src/models/gcc_solver.cpp
    src/models/gcc_solver_pcg.cpp
    src/models/gcc_preconditioner.cpp
    ${VCF_SOURCES}
)
target_link_libraries(test_cross_operator PRIVATE CosmicCommon OpenMP::OpenMP_CXX)
target_include_directories(test_cross_operator PRIVATE
    src
    src/models
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)
add_test(NAME cross_operator_unit COMMAND test_cross_operator)
set_tests_properties(cross_operator_unit PROPERTIES
    LABELS "unit;ssgwas;apy;e3d2d" TIMEOUT 300
    ENVIRONMENT "OMP_NUM_THREADS=1;OPENBLAS_NUM_THREADS=1;MKL_NUM_THREADS=1"
)

add_executable(test_apy_ginv_d15 tests/test_apy_ginv_d15.cpp)
target_link_libraries(test_apy_ginv_d15 PRIVATE CosmicCommon)
add_test(NAME apy_ginv_eq8_schur_equivalence COMMAND test_apy_ginv_d15)
set_tests_properties(apy_ginv_eq8_schur_equivalence PROPERTIES
    LABELS "unit;ssgwas;apy;d15" TIMEOUT 120)

add_executable(eval_repeatability_reml tools/eval_repeatability_reml.cpp)
target_link_libraries(eval_repeatability_reml PRIVATE CosmicCommon)
target_include_directories(eval_repeatability_reml PRIVATE
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)
endif()

install(TARGETS cosmicscanner RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
