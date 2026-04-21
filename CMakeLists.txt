cmake_minimum_required(VERSION 3.15)
project(bardmux VERSION 2.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(
    -O3
    -march=native
    -funroll-loops
    -fomit-frame-pointer
    -fPIC
    $<$<CONFIG:Release>:-DNDEBUG>
)

# Link-time optimisation for Release builds
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# Only zlib and pthreads are required — kseq.h is vendored in src/
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)

set(SOURCES
    src/main.cpp
    src/fastq_reader.cpp
    src/anchor_finder.cpp
    src/barcode_matcher.cpp
    src/report_writer.cpp
    src/thread_pool.cpp
)

add_executable(bardmux ${SOURCES})

target_include_directories(bardmux PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)

target_link_libraries(bardmux PRIVATE
    ZLIB::ZLIB
    Threads::Threads
)

enable_testing()

install(TARGETS bardmux DESTINATION bin)

message(STATUS "bardmux v2.0.0  LTO=${CMAKE_INTERPROCEDURAL_OPTIMIZATION}")
