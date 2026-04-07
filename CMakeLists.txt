cmake_minimum_required(VERSION 3.15)
project(CosmicScanner CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB SCANNER_SOURCES src/*.cpp)

add_executable(cosmicscanner
  ${SCANNER_SOURCES}
)

target_include_directories(cosmicscanner PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_link_libraries(cosmicscanner PRIVATE CosmicCommon)

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

install(TARGETS cosmicscanner RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
