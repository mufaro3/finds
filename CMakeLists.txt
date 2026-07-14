cmake_minimum_required(VERSION 3.20)

project(
    Simulation
    VERSION 1.0
    LANGUAGES C
)

# Compiler settings

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Dependencies

find_package(OpenCL REQUIRED)
find_package(GSL REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(CHECK REQUIRED check)

# Third-party libraries

# Header-only gnuplot_i library
add_library(gnuplot_i INTERFACE)

target_include_directories(gnuplot_i
    INTERFACE
        ${CMAKE_CURRENT_SOURCE_DIR}/lib
)

# Simulation library

file(GLOB_RECURSE LIB_SOURCES
    src/*.c
)

add_library(simulation STATIC
    ${LIB_SOURCES}
)

target_include_directories(simulation
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(simulation
    PUBLIC
        OpenCL::OpenCL
        GSL::gsl
        GSL::gslcblas
        gnuplot_i
        m
)

target_compile_definitions(simulation
    PUBLIC
        CL_TARGET_OPENCL_VERSION=300
)

# Programs

file(GLOB PROGRAM_SOURCES
    programs/*.c
)

foreach(PROGRAM ${PROGRAM_SOURCES})

    get_filename_component(PROGRAM_NAME ${PROGRAM} NAME_WE)

    add_executable(${PROGRAM_NAME}
        ${PROGRAM}
    )

    target_link_libraries(${PROGRAM_NAME}
        PRIVATE
            simulation
    )

endforeach()

# Tests

enable_testing()

file(GLOB TEST_SOURCES
    tests/*.c
)

add_executable(unit_tests
    ${TEST_SOURCES}
)

target_include_directories(unit_tests
    PRIVATE
        ${CHECK_INCLUDE_DIRS}
)

target_link_libraries(unit_tests
    PRIVATE
        simulation
        ${CHECK_LIBRARIES}
)

add_test(
    NAME unit_tests
    COMMAND unit_tests
)

# Documentation

find_package(Doxygen)

if(DOXYGEN_FOUND)

    set(DOXYGEN_GENERATE_HTML YES)
    set(DOXYGEN_GENERATE_LATEX NO)

    doxygen_add_docs(
        docs
        include
        src
    )

endif()
