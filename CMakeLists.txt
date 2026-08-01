cmake_minimum_required(VERSION 3.15)
project(MCMR LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(pybind11 REQUIRED)

set(MCMR_SOURCES
    src/main.cpp
    src/geometry.cpp
    src/physics.cpp
    src/simulation.cpp
    src/io.cpp
)

pybind11_add_module(_mcmr_cpp ${MCMR_SOURCES})
target_include_directories(_mcmr_cpp PRIVATE include)

install(TARGETS _mcmr_cpp DESTINATION mcmr)

install(DIRECTORY mcmr/ DESTINATION mcmr)

file(GLOB HTML_FILES "${CMAKE_CURRENT_SOURCE_DIR}/data/*.html")
install(FILES ${HTML_FILES} DESTINATION data)