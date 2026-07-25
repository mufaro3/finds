cmake_minimum_required(VERSION 3.20)
project(FINDS VERSION 1.0 LANGUAGES C)

# ------------------------------------------------------------------------------
# Compiler settings
# ------------------------------------------------------------------------------

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# ------------------------------------------------------------------------------
# Dependencies
# ------------------------------------------------------------------------------

find_package(OpenCL REQUIRED)
find_package(GSL REQUIRED)

find_package(PkgConfig REQUIRED)
pkg_check_modules(CHECK REQUIRED check)
find_package(Curses REQUIRED)

find_package(HDF5 REQUIRED COMPONENTS C)

# ------------------------------------------------------------------------------
# Third-party libraries
# ------------------------------------------------------------------------------
set(LIB_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# gnuplot_i
set(GNUPLOT_ROOT ${LIB_DIR}/gnuplot_i)
add_library(gnuplot_i STATIC ${GNUPLOT_ROOT}/gnuplot_i.c)
target_include_directories(gnuplot_i PUBLIC ${LIB_DIR})
target_link_libraries(gnuplot_i PUBLIC m)

# progressbar
set(PROGRESSBAR_ROOT ${LIB_DIR}/progressbar)
add_subdirectory(${PROGRESSBAR_ROOT})
target_include_directories(progressbar
    PRIVATE
    ${PROGRESSBAR_ROOT}/include

    PUBLIC
    ${CURSES_INCLUDE_DIRS}
)

target_link_libraries(progressbar PUBLIC ${CURSES_LIBRARIES})

# tomlc17
set(TOMLC17_ROOT ${LIB_DIR}/tomlc17)
add_library(tomlc17 STATIC ${TOMLC17_ROOT}/src/tomlc17.c)
target_include_directories(tomlc17 PUBLIC ${TOMLC17_ROOT}/src)

# ------------------------------------------------------------------------------
# Simulation library
# ------------------------------------------------------------------------------

file(GLOB_RECURSE LIB_SOURCES CONFIGURE_DEPENDS src/*.c)

add_library(finds STATIC ${LIB_SOURCES})

target_include_directories(finds
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include

    PRIVATE
    ${OpenCL_INCLUDE_DIRS}
    ${GSL_INCLUDE_DIRS}
    ${CHECK_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${PROGRESSBAR_ROOT}/include
)

target_link_libraries(finds PUBLIC
    OpenCL::OpenCL
    GSL::gsl
    GSL::gslcblas
    gnuplot_i
    progressbar
    tomlc17
    HDF5::HDF5
    m
)

target_compile_definitions(finds PUBLIC CL_TARGET_OPENCL_VERSION=300)
target_compile_options(finds PRIVATE -Wall -Wextra -Wno-missing-braces -O3)

# ------------------------------------------------------------------------------
# Programs
# ------------------------------------------------------------------------------

file(GLOB PROGRAM_SOURCES CONFIGURE_DEPENDS programs/*.c)

foreach(PROGRAM ${PROGRAM_SOURCES})
    get_filename_component(PROGRAM_NAME ${PROGRAM} NAME_WE)
    add_executable(${PROGRAM_NAME} ${PROGRAM})
    target_include_directories(${PROGRAM_NAME} PRIVATE
        ${PROGRESSBAR_ROOT}/include
        ${TOMLC17_ROOT}/src
    )
    target_link_libraries(${PROGRAM_NAME} PRIVATE finds)
    target_compile_options(${PROGRAM_NAME} PRIVATE -Wall -Wextra
        -Wno-missing-braces -lm -O3)
endforeach()

# ------------------------------------------------------------------------------
# Tests
# ------------------------------------------------------------------------------

enable_testing()

file(GLOB TEST_SOURCES CONFIGURE_DEPENDS tests/*.c)

add_executable(unit_tests ${TEST_SOURCES})

target_include_directories(unit_tests PRIVATE ${CHECK_INCLUDE_DIRS})
target_link_libraries(unit_tests PRIVATE finds ${CHECK_LIBRARIES})

add_test(NAME unit_tests COMMAND unit_tests)

# ------------------------------------------------------------------------------
# Documentation
# ------------------------------------------------------------------------------

find_package(Doxygen)
find_package(LATEX)

if(DOXYGEN_FOUND)

    set(DOXYGEN_GENERATE_HTML YES)
    set(DOXYGEN_GENERATE_LATEX YES)

    doxygen_add_docs(docs include src)

    if(LATEX_FOUND)
        add_custom_command(
            TARGET docs POST_BUILD
            COMMAND ${CMAKE_MAKE_PROGRAM}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/latex
            COMMENT "Building Doxygen PDF documentation"
        )
    endif()

endif()
