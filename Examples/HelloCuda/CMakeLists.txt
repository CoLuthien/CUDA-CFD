# Target name
set(target HelloCuda)

# Includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR})


# Sources
file(GLOB sources
    ${CMAKE_CURRENT_SOURCE_DIR}/*.f90)

# Build executable
add_executable(${target}
    ${sources})
target_include_directories(${target}
    PRIVATE
    ${CMAKE_Fortran_MODULE_DIRECTORY}
)

# Project options
#set_target_properties(${target}
#    PROPERTIES
#    ${DEFAULT_PROJECT_OPTIONS}
#)

# Compile options
if (CMAKE_CXX_COMPILER_ID MATCHES "MSVC")
    add_definitions(-D_USE_MATH_DEFINES) # for M_PI
endif ()

target_compile_options(${target}
    PRIVATE

    PUBLIC
    ${DEFAULT_COMPILE_OPTIONS}
    INTERFACE
)


# Link libraries
target_link_libraries(${target}
    PUBLIC
    FortFlow
    PRIVATE
    ${DEFAULT_LINKER_OPTIONS}
    )

message(${CMAKE_Fortran_COMPILER_ID})