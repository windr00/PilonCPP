cmake_minimum_required(VERSION 3.10)
project(champhash VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Header-only interface library
add_library(champhash INTERFACE)
add_library(champhash::champhash ALIAS champhash)

target_include_directories(champhash INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# Tests
option(BUILD_TESTS "Build champhash tests" ON)
if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
