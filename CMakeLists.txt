cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(nuanalysis)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

set(CMAKE_CXX_STANDARD 20)
find_package(nlohmann_json 3.2.0 REQUIRED)
find_package(ROOT REQUIRED)
include("${ROOT_USE_FILE}")
include_directories(
    ${ROOT_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_library(event STATIC src/event.cxx)
# add_library(analysis STATIC src/analysis.cxx)
# add_library(analysis_plot STATIC src/analysis_plot.cxx)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -Wextra -Wpedantic")
add_executable(nuanalysis main.cxx)
#target_compile_options(nuanalysis PUBLIC -pg)
#target_link_options(nuanalysis PUBLIC -pg)
target_link_libraries(nuanalysis ROOT::Tree ROOT::EG event nlohmann_json::nlohmann_json)
install(TARGETS nuanalysis DESTINATION bin)

add_executable(pdkanalysis pdkanalysis.cxx)
target_link_libraries(pdkanalysis ROOT::Tree ROOT::EG event nlohmann_json::nlohmann_json)
install(TARGETS pdkanalysis DESTINATION bin)
