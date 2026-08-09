cmake_minimum_required(VERSION 3.14)
project(piloncpp VERSION 1.2.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ---- Static htslib (downloaded and built on first configure) ----
include(ExternalProject)

set(HTSLIB_SRC ${CMAKE_BINARY_DIR}/htslib_src)

# htslib's autoconf defaults to CFLAGS="-g -Wall -O2"; override to drop -g so
# debug info from the static lib does not leak into the final binary.
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(HTSLIB_CFLAGS "-O2 -DNDEBUG")
else()
    set(HTSLIB_CFLAGS "-g -O2")
endif()

ExternalProject_Add(htslib_ep
    URL https://github.com/samtools/htslib/releases/download/1.19/htslib-1.19.tar.bz2
    SOURCE_DIR ${HTSLIB_SRC}
    BUILD_IN_SOURCE ON
    CONFIGURE_COMMAND env CFLAGS=${HTSLIB_CFLAGS} ${HTSLIB_SRC}/configure --disable-libcurl --disable-plugins --without-libdeflate --disable-bz2 --disable-lzma
    BUILD_COMMAND make -C ${HTSLIB_SRC} lib-static
    INSTALL_COMMAND ""
    BUILD_BYPRODUCTS ${HTSLIB_SRC}/libhts.a
)

# Create an imported target for the static library
add_library(hts STATIC IMPORTED)
set_target_properties(hts PROPERTIES IMPORTED_LOCATION ${HTSLIB_SRC}/libhts.a)
add_dependencies(hts htslib_ep)

# ---- Include paths ----
include_directories(${CMAKE_SOURCE_DIR}/include)
include_directories(${HTSLIB_SRC}/htslib)

# ---- Source files ----
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
    src/scaffold.cpp
    src/main.cpp
)

add_executable(piloncpp ${SOURCES})

# ---- Link (static htslib + static C++ runtime) ----
find_package(Threads REQUIRED)

target_link_libraries(piloncpp PRIVATE
    hts
    Threads::Threads
    dl m z
)

target_link_options(piloncpp PRIVATE
    -static-libstdc++
    -static-libgcc
)

# ---- Compiler flags ----
target_compile_options(piloncpp PRIVATE -O2 -Wall -Wextra)

# ---- Install ----
install(TARGETS piloncpp DESTINATION bin)
