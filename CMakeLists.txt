# CMake version
cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
cmake_policy(SET CMP0063 NEW) # visibility

# Include cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/Builds/CMake")

# Declare project
project(FortFlow)

# Useful paths
set(LIBRARY_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty)
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/build/mod)
enable_language(Fortran)

# Set output directories
set(DEFAULT_CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)


# Set enable output of compile commands during generation
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_Fortran_FLAGS
    "-fopenmp -check bound")



# Build type - Release by default
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
elseif (CMAKE_BUILD_TYPE MATCHES "Debug")
    add_definitions(-DDEBUG)
endif()
message(STATUS "CMake build type: ${CMAKE_BUILD_TYPE}")

# CUDA
option(USE_CUDA "Use CUDA features" ON)
set(CUDA_CRT_LINKAGE "dynamic" CACHE STRING
                     "CUDA targets' CRT options. For 'static', use /MT or /MTd. For 'dynamic', use /MD or /MDd")
if (USE_CUDA)
    set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)
    set(CUDA_SEPARABLE_COMPILATION OFF)
    set(CUDA_STANDARD 17)
    find_package(CUDA)
    if (CUDA_FOUND) # Modify project options for CUDA build
        enable_language(CUDA)
        cuda_select_nvcc_arch_flags(CUDA_NVCC_FLAGS "7.5" "Turing")
        # We will use CUDA_TARGET_OPTIONS instead of CUDA_NVCC_FLAGS to be explicit.
        # start with architecture options
        list(APPEND CUDA_TARGET_OPTIONS
            "--std=c++17" "-Wno-deprecated-gpu-targets")
        if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
            list(APPEND CUDA_TARGET_OPTIONS "-Xcompiler" "/FS") # prevent race on PDB
            if (CUDA_CRT_LINKAGE STREQUAL "dynamic")
                if (CMAKE_BUILD_TYPE MATCHES Debug)
                    list(APPEND CUDA_TARGET_OPTIONS "-Xcompiler" "/MDd" "--debug" "-O0")
                else()
                    list(APPEND CUDA_TARGET_OPTIONS "-Xcompiler" "/MD" "-O2")
                endif()
            elseif (CUDA_CRT_LINKAGE STREQUAL "static")
                if(CMAKE_BUILD_TYPE MATCHES Debug)
                    list(APPEND CUDA_TARGET_OPTIONS "-Xcompiler" "/MTd" "--debug" "-O0")
                else()
                    list(APPEND CUDA_TARGET_OPTIONS "-Xcompiler" "/MT" "-O2")
                endif()
            endif()
            list(APPEND DEFAULT_COMPILE_OPTIONS
                /wd4515 # -> disable warning: [Thrust issue] 'detail': namespace uses itself
                /wd4996 # -> disable warning: [Thrust issue] 'strerror': This function or variable may be unsafe.
                /wd4800 # -> disable warning: [Thrust issue] forcing value to bool 'true' or 'false' (performance warning)
            )
        endif()
        add_compile_definitions(CUBBYFLOW_USE_CUDA)
        message(STATUS "Using CUDA: ${CUDA_VERSION}")
        message(STATUS "  CRT: ${CUDA_CRT_LINKAGE}")
        message(STATUS "  options: ${CUDA_TARGET_OPTIONS}")
        message(STATUS "  nvcc flags: ${CUDA_NVCC_FLAGS}")
        # note: This is for CI build. It fails to find this library so we are mentioning it here
        if(EXISTS ${CUDA_cudadevrt_LIBRARY})
            message(STATUS "  cudadevrt: ${CUDA_cudadevrt_LIBRARY}")
            get_filename_component(CUDA_LIB_DIR ${CUDA_cudadevrt_LIBRARY} DIRECTORY)
            link_directories(${CUDA_LIB_DIR})
            message(STATUS "  lib: ${CUDA_LIB_DIR}")
        endif()
    else()
        set(USE_CUDA OFF)
    endif()
endif()

# Find TBB
#include(Builds/CMake/FindTBB.cmake)

# Tasking system options
#include(Builds/CMake/TaskingSystemOptions.cmake)

# Compile options
#include(Builds/CMake/CompileOptions.cmake)

# Code coverage - Debug only
# NOTE: Code coverage results with an optimized (non-Debug) build may be misleading
option(BUILD_COVERAGE "Build code coverage" OFF)
if (CMAKE_BUILD_TYPE MATCHES Debug AND CMAKE_COMPILER_IS_GNUCXX AND BUILD_COVERAGE)
    include(CodeCoverage)
    setup_target_for_coverage(${PROJECT_NAME}_coverage UnitTests coverage)
endif()



# Project modules

## ThirdParty modules
add_subdirectory(Sources/Core)

add_subdirectory(Examples/HelloCuda)
if (USE_CUDA)
    #add_subdirectory(Examples/HelloCuda)
endif()

#if (NOT BUILD_SONARCLOUD)
#    add_subdirectory(Tests/MemPerfTests)
#    add_subdirectory(Tests/TimePerfTests)
#    add_subdirectory(Examples/HelloFluidSim)
#    add_subdirectory(Examples/HybridLiquidSim)
#    add_subdirectory(Examples/LevelSetLiquidSim)
#    add_subdirectory(Examples/Obj2Sdf)
#    add_subdirectory(Examples/Particles2Obj)
#    add_subdirectory(Examples/Particles2Xml)
#    add_subdirectory(Examples/SmokeSim)
#    add_subdirectory(Examples/SPHSim)
#endif()
