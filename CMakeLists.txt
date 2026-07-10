cmake_minimum_required(VERSION 3.20)

project(Simulation
    VERSION 1.0
    LANGUAGES C
)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(OpenCL REQUIRED)
find_package(GSL REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(CHECK REQUIRED check)

include_directories(
    include
    ${OpenCL_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS}
    ${CHECK_INCLUDE_DIRS}
)

file(GLOB_RECURSE SOURCES
    src/*.c
)

add_executable(simulation
    ${SOURCES}
)

target_link_libraries(simulation
    OpenCL::OpenCL
    GSL::gsl
    GSL::gslcblas
)

enable_testing()

file(GLOB TEST_SOURCES tests/*.c)

add_executable(unit_tests
    ${TEST_SOURCES}
)

target_link_libraries(unit_tests
    OpenCL::OpenCL
    GSL::gsl
    GSL::gslcblas
    ${CHECK_LIBRARIES}
)

target_compile_definitions(simulation PRIVATE
    CL_TARGET_OPENCL_VERSION=300
)

target_compile_definitions(unit_tests PRIVATE
    CL_TARGET_OPENCL_VERSION=300
)

add_test(
    NAME unit_tests
    COMMAND unit_tests
)

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
