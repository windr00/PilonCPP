cmake_minimum_required(VERSION 3.14)
project(piloncpp VERSION 1.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find htslib directly (Ubuntu libhts-dev doesn't provide pkg-config .pc file)
find_path(HTSLIB_INCLUDE_DIR hts.h PATH_SUFFIXES htslib)
find_library(HTSLIB_LIBRARY hts)

if(NOT HTSLIB_INCLUDE_DIR OR NOT HTSLIB_LIBRARY)
    message(FATAL_ERROR "htslib not found. Please install libhts-dev.")
endif()

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${HTSLIB_INCLUDE_DIR})

# Source files
set(SOURCES
    src/bases.cpp
    src/basesum.cpp
    src/pileup.cpp
    src/pilon.cpp
    src/genome.cpp
    src/vcf.cpp
    src/assembler.cpp
    src/gapfiller.cpp
    src/normal.cpp
    src/bamfile.cpp
    src/utils.cpp
    src/main.cpp
)

# Executable
add_executable(piloncpp ${SOURCES})

# Link libraries
find_library(BZ2_LIBRARIES bz2)
find_library(Z_LIBRARIES z)
find_library(CURL_LIBRARIES curl)
find_package(Threads REQUIRED)

target_link_libraries(piloncpp ${HTSLIB_LIBRARY} Threads::Threads)
if(BZ2_LIBRARIES)
    target_link_libraries(piloncpp ${BZ2_LIBRARIES})
endif()
if(Z_LIBRARIES)
    target_link_libraries(piloncpp ${Z_LIBRARIES})
endif()
if(CURL_LIBRARIES)
    target_link_libraries(piloncpp ${CURL_LIBRARIES})
endif()

# Compiler options
target_compile_options(piloncpp PRIVATE -O2 -Wall -Wextra)

# Install
install(TARGETS piloncpp DESTINATION bin)
