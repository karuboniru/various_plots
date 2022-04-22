cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(catch_abnormal_wtrue)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

set(CMAKE_CXX_STANDARD 20)
find_package(ROOT REQUIRED)
include("${ROOT_USE_FILE}")
include_directories(
    ${ROOT_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_library(event STATIC src/event.cxx)
# add_library(analysis STATIC src/analysis.cxx)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -Wextra -Wpedantic")
add_executable(pion main.cxx)
target_link_libraries(pion ROOT::Tree ROOT::EG event)
install(TARGETS pion DESTINATION bin)

add_executable(compare compare.cxx)
target_link_libraries(compare ROOT::Tree ROOT::EG event)
install(TARGETS compare DESTINATION bin)

# add_executable(compare_cut compare_cut.cxx)
# target_link_libraries(compare_cut ROOT::Tree ROOT::EG event)
# install(TARGETS compare_cut DESTINATION bin)

# add_executable(muon muon.cxx)
# target_link_libraries(muon ROOT::Tree ROOT::EG event )
# install(TARGETS muon DESTINATION bin)