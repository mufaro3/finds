find_package(pybind11 REQUIRED)
pybind11_add_module(finds
    pyfinds.cpp
    system.cpp
    datastream.cpp
    simulation.cpp
)
target_link_libraries(finds PRIVATE finds_core)
target_include_directories(finds PRIVATE ../include)
