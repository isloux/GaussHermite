cmake_minimum_required(VERSION 3.20)

project(GaussHermite
    VERSION 1.0
    DESCRIPTION "Fast Gauss-Hermite Quadrature Implementation"
    LANGUAGES CXX)

# Set C++20 standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Compiler warnings and optimizations
if(MSVC)
    add_compile_options(/W4 /WX)
else()
    add_compile_options(
        -Wall -Wextra -Wpedantic -Werror
        -Wconversion -Wsign-conversion
        -Wcast-qual -Wcast-align
        -Wold-style-cast
        -O3 -march=native
    )
endif()

# Main executable
add_executable(gauss_hermite gauss_hermite.cpp)

# Link against math library on Unix
if(UNIX AND NOT APPLE)
    target_link_libraries(gauss_hermite PRIVATE m)
endif()

# Optional: Enable link-time optimization for release builds
include(CheckIPOSupported)
check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
if(ipo_supported)
    set_property(TARGET gauss_hermite PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
