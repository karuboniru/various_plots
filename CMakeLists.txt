cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(pions)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

set(CMAKE_CXX_STANDARD 20)
find_package(ROOT REQUIRED)
include("${ROOT_USE_FILE}")
include_directories(
    ${ROOT_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_library(event STATIC src/event.cxx)
add_library(analysis STATIC src/analysis.cxx)
add_library(analysis_plot STATIC src/analysis_plot.cxx)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -Wextra -Wpedantic")
add_executable(pion main.cxx)
target_link_libraries(pion ROOT::Tree ROOT::EG analysis analysis_plot event)
install(TARGETS pion DESTINATION bin)

add_executable(compare compare.cxx)
target_link_libraries(compare ROOT::Tree ROOT::EG analysis analysis_plot event)
install(TARGETS compare DESTINATION bin)
