cmake_minimum_required(VERSION 3.16)

project(cpm_module LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(ROOT REQUIRED COMPONENTS Core Hist Physics RIO Tree)

add_library(cpm_core
  CPMAverageCorrectionReconstruction.cc
  CPMReconstructionHelper.cc
  CPMVoxelContainerv1.cc
)

target_include_directories(cpm_core
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(cpm_core
  PUBLIC
    ROOT::Core
    ROOT::Hist
    ROOT::Physics
    ROOT::RIO
    ROOT::Tree
)

add_executable(test_local_line_poca
  test_local_line_poca.cc
)

target_link_libraries(test_local_line_poca
  PRIVATE cpm_core
)

enable_testing()
add_test(NAME local_line_poca COMMAND test_local_line_poca)
