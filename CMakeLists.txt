cmake_minimum_required(VERSION 3.15)
project(CosmicScanner VERSION 1.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

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

add_executable(test_sparse_reml tests/test_sparse_reml.cpp)
target_link_libraries(test_sparse_reml PRIVATE CosmicCommon)

add_executable(test_multi_vc tests/test_multi_vc.cpp)
target_link_libraries(test_multi_vc PRIVATE CosmicCommon)

add_executable(test_spa_accuracy tests/test_spa_accuracy.cpp src/spa_test.cpp)
target_link_libraries(test_spa_accuracy CosmicCommon)
target_include_directories(test_spa_accuracy PRIVATE src)

add_executable(test_regenie tests/test_regenie.cpp src/ridge_regenie.cpp)
target_link_libraries(test_regenie CosmicCommon)
target_include_directories(test_regenie PRIVATE src)

add_executable(test_input_adaptation_units
    tests/test_input_adaptation_units.cpp
    src/input/input_common.cpp
    src/input/input_source_reader.cpp
    src/input/time_normalizer.cpp
    src/input/wide_to_long.cpp
)
target_include_directories(test_input_adaptation_units PRIVATE src)

add_executable(eval_repeatability_reml tools/eval_repeatability_reml.cpp)
target_link_libraries(eval_repeatability_reml PRIVATE CosmicCommon)
target_include_directories(eval_repeatability_reml PRIVATE
    ${CMAKE_SOURCE_DIR}/libCosmicCore/src
)

install(TARGETS cosmicscanner RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
