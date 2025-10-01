cmake_minimum_required(VERSION 3.12)
project(phmm-tree VERSION 1.0 LANGUAGES C CXX)

# Set C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find OpenMP
find_package(OpenMP)

# Source files
set(SOURCES
    fitch.c
    phylip.c
    dist.c
    neighbor.c
    upgma.c
    kitsch.c
    class_functions.cpp
    process_matrices.cpp
    public_functions.cpp
    process_alignments.cpp
    process_hmms.cpp
    process_usearch.cpp
    process_prc.cpp
    HMMTree.cpp
    phylip_draw_tree.cpp
    hhsuite.cpp
)

# Create executable
add_executable(phmm-tree ${SOURCES})

# Compiler-specific optimizations
if(CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID STREQUAL "Clang")
    target_compile_options(phmm-tree PRIVATE -O3 -Wall)
endif()

# Link math library
target_link_libraries(phmm-tree PRIVATE m)

# Add OpenMP support if available
if(OpenMP_C_FOUND)
    message(STATUS "OpenMP found - enabling parallel compilation")
    target_link_libraries(phmm-tree PRIVATE OpenMP::OpenMP_C)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(phmm-tree PRIVATE OpenMP::OpenMP_CXX)
    endif()
    target_compile_definitions(phmm-tree PRIVATE OPENMP_ENABLED)
else()
    message(WARNING "OpenMP not found - compiling sequential version")
endif()

# Add thread detection flag
target_compile_definitions(phmm-tree PRIVATE AUTO_THREAD_DETECTION)

# Installation
install(TARGETS phmm-tree DESTINATION bin)

# Print configuration
message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "C++ Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
if(OpenMP_C_FOUND)
    message(STATUS "OpenMP version: ${OpenMP_C_VERSION}")
endif()
