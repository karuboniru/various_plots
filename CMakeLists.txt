cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(nuanalysis)
set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)

set(CMAKE_CXX_STANDARD 20)
# add_compile_options(-fconcepts)
find_package(nlohmann_json 3.2.0 REQUIRED)
find_package(ROOT REQUIRED)
include("${ROOT_USE_FILE}")
include_directories(
    ${ROOT_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

add_library(event STATIC src/event.cxx)
add_library(tools STATIC src/tools.cxx)
# add_library(analysis STATIC src/analysis.cxx)
# add_library(analysis_plot STATIC src/analysis_plot.cxx)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wall -Wextra -Wpedantic -flto -fconcepts")

add_executable(rdf-rho rdf-rho.cxx)
target_link_libraries(rdf-rho ROOT::Tree ROOT::EG ROOT::ROOTDataFrame event tools nlohmann_json::nlohmann_json)
install(TARGETS rdf-rho DESTINATION bin)

add_executable(rdf-eta rdf-eta.cxx)
target_link_libraries(rdf-eta ROOT::Tree ROOT::EG ROOT::ROOTDataFrame event tools nlohmann_json::nlohmann_json)
install(TARGETS rdf-eta DESTINATION bin)

add_executable(rdf-pion rdf-pion.cxx)
target_link_libraries(rdf-pion ROOT::Tree ROOT::EG ROOT::ROOTDataFrame event tools nlohmann_json::nlohmann_json)
install(TARGETS rdf-pion DESTINATION bin)

add_executable(rdf-pdp rdf-pdp.cxx)
target_link_libraries(rdf-pdp ROOT::Tree ROOT::EG ROOT::ROOTDataFrame event tools nlohmann_json::nlohmann_json)
install(TARGETS rdf-pdp DESTINATION bin)

add_executable(rdf-omega rdf-omega.cxx)
target_link_libraries(rdf-omega ROOT::Tree ROOT::EG ROOT::ROOTDataFrame event tools nlohmann_json::nlohmann_json)
install(TARGETS rdf-omega DESTINATION bin)

add_executable(genie_kaon genie_kaon.cxx)
target_link_libraries(genie_kaon ROOT::Tree ROOT::EG ROOT::ROOTDataFrame event tools nlohmann_json::nlohmann_json)
install(TARGETS genie_kaon DESTINATION bin)

add_executable(channel_QW_plot channel_QW_plot.cxx)
target_link_libraries(channel_QW_plot ROOT::Tree ROOT::EG ROOT::ROOTDataFrame event tools nlohmann_json::nlohmann_json)
install(TARGETS channel_QW_plot DESTINATION bin)

add_executable(minerva_pi_plus minerva_pi_plus.cxx)
target_link_libraries(minerva_pi_plus ROOT::Tree ROOT::EG ROOT::ROOTDataFrame event tools nlohmann_json::nlohmann_json)
install(TARGETS minerva_pi_plus DESTINATION bin)

add_executable(minerva_pi_plus_nofsi minerva_pi_plus_nofsi.cxx)
target_link_libraries(minerva_pi_plus_nofsi ROOT::Tree ROOT::EG ROOT::ROOTDataFrame event tools nlohmann_json::nlohmann_json)
install(TARGETS minerva_pi_plus_nofsi DESTINATION bin)

add_executable(minerva_pi_plus_nofsi_1pi minerva_pi_plus_nofsi_1pi.cxx)
target_link_libraries(minerva_pi_plus_nofsi_1pi ROOT::Tree ROOT::EG ROOT::ROOTDataFrame event tools nlohmann_json::nlohmann_json)
install(TARGETS minerva_pi_plus_nofsi_1pi DESTINATION bin)

add_executable(nuwro_plot nuwro_plot.cxx)
target_link_libraries(nuwro_plot ROOT::Tree ROOT::EG ROOT::ROOTDataFrame event tools nlohmann_json::nlohmann_json)
install(TARGETS nuwro_plot DESTINATION bin)

add_executable(stack_plot stack_plot.cxx)
target_link_libraries(stack_plot ROOT::Tree ROOT::EG ROOT::ROOTDataFrame event tools nlohmann_json::nlohmann_json)
install(TARGETS stack_plot DESTINATION bin)

add_executable(GENIE_all_compare GENIE_all_compare.cxx)
target_link_libraries(GENIE_all_compare ROOT::Tree ROOT::EG ROOT::ROOTDataFrame event tools nlohmann_json::nlohmann_json)
install(TARGETS GENIE_all_compare DESTINATION bin)


add_executable(genie_proton genie_proton.cxx)
target_link_libraries(genie_proton ROOT::Tree ROOT::EG ROOT::ROOTDataFrame event tools nlohmann_json::nlohmann_json)
install(TARGETS genie_proton DESTINATION bin)

add_executable(trace_particle trace_particle.cxx)
target_link_libraries(trace_particle ROOT::Tree ROOT::EG ROOT::ROOTDataFrame event tools nlohmann_json::nlohmann_json)
install(TARGETS trace_particle DESTINATION bin)