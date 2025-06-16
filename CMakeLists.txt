cmake_minimum_required(VERSION 3.10)
project(bellman_ford_gaspi C)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Use MPI compiler wrapper (mpicc)
set(CMAKE_C_COMPILER mpicc)

# Get GASPI environment variables
if(NOT DEFINED ENV{GPI_INC_DIR} OR NOT DEFINED ENV{GPI_LIB_DIR})
    message(WARNING "GASPI environment variables not set. Trying to source...")
    execute_process(COMMAND bash -c "source /opt/nfs/config/source_gaspi.sh > /dev/null 2>&1 && env"
                    OUTPUT_VARIABLE GASPI_ENV)
    string(REGEX MATCH "GPI_INC_DIR=([^\n]+)" GPI_INC_DIR_MATCH ${GASPI_ENV})
    string(REGEX MATCH "GPI_LIB_DIR=([^\n]+)" GPI_LIB_DIR_MATCH ${GASPI_ENV})
    
    if(GPI_INC_DIR_MATCH)
        set(GPI_INC_DIR ${CMAKE_MATCH_1})
    else()
        message(FATAL_ERROR "Could not find GPI_INC_DIR. Please source GASPI environment first.")
    endif()
    
    if(GPI_LIB_DIR_MATCH)
        set(GPI_LIB_DIR ${CMAKE_MATCH_1})
    else()
        message(FATAL_ERROR "Could not find GPI_LIB_DIR. Please source GASPI environment first.")
    endif()
else()
    set(GPI_INC_DIR $ENV{GPI_INC_DIR})
    set(GPI_LIB_DIR $ENV{GPI_LIB_DIR})
endif()

message(STATUS "GASPI include directory: ${GPI_INC_DIR}")
message(STATUS "GASPI library directory: ${GPI_LIB_DIR}")

# Verify gaspi.h exists
if(NOT EXISTS "${GPI_INC_DIR}/GASPI.h")
    message(FATAL_ERROR "gaspi.h not found in ${GPI_INC_DIR}!")
endif()

# Verify library exists
if(NOT EXISTS "${GPI_LIB_DIR}/libGPI2.so")
    message(FATAL_ERROR "libGPI2.so not found in ${GPI_LIB_DIR}!")
endif()

add_executable(bellman_ford_gaspi
    src/bellman-ford-gaspi.c
    src/graph.c
)

# Include directories: 
# - Current source directory (.)
# - Project include directory
# - GASPI include directory
target_include_directories(bellman_ford_gaspi PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_SOURCE_DIR}/include"
    "${GPI_INC_DIR}"
    ${CUDA_INCLUDE_DIRS}
)

# Find the CUDA driver library
find_library(CUDA_DRIVER_LIBRARY
  NAMES cuda
  HINTS ${CUDA_TOOLKIT_ROOT_DIR}/lib64
)

if (NOT CUDA_DRIVER_LIBRARY)
  message(FATAL_ERROR "Could not find libcuda.so! You need the CUDA driver library on your link path.")
endif()

# Link directories and libraries
target_link_directories(bellman_ford_gaspi PRIVATE
    "${GPI_LIB_DIR}"
)

target_link_libraries(bellman_ford_gaspi PRIVATE 
    GPI2 pthread m
    ${CUDA_DRIVER_LIBRARY}
)

# Additional compiler flags
target_compile_options(bellman_ford_gaspi PRIVATE
    -Wall -Wextra
)

# Add compilation definitions
target_compile_definitions(bellman_ford_gaspi PRIVATE
    -D_GNU_SOURCE
    -DGASPI_EXT_PASSIVE=1
)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set(CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -Wl,-rpath=${CUDA_TOOLKIT_ROOT_DIR}/lib64"
)