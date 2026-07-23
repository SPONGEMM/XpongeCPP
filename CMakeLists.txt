cmake_minimum_required(VERSION 3.18)
project(XpongeCPP LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(PYBIND11_FINDPYTHON ON)

find_package(pybind11 CONFIG REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(HighFive CONFIG QUIET)
if(NOT HighFive_FOUND)
    include(FetchContent)
    set(HIGHFIVE_FIND_HDF5 OFF CACHE BOOL "Do not let vendored HighFive find HDF5" FORCE)
    set(HIGHFIVE_UNIT_TESTS OFF CACHE BOOL "Disable vendored HighFive tests" FORCE)
    set(HIGHFIVE_EXAMPLES OFF CACHE BOOL "Disable vendored HighFive examples" FORCE)
    set(HIGHFIVE_BUILD_DOCS OFF CACHE BOOL "Disable vendored HighFive documentation" FORCE)
    FetchContent_Declare(
        highfive
        URL https://github.com/highfive-devs/highfive/archive/refs/tags/v3.3.0.tar.gz
        URL_HASH SHA256=325cfbcf0c0296a6dd26f3b088801b7ebb8d6f109c0565c11d2d8c4af3253bff
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_MakeAvailable(highfive)
    set_property(DIRECTORY "${highfive_SOURCE_DIR}" PROPERTY EXCLUDE_FROM_ALL TRUE)
endif()

pybind11_add_module(_core
    cpp/assign/assign.cpp
    cpp/assign/bond_order.cpp
    cpp/assign/gaff.cpp
    cpp/assign/mol2_assignment.cpp
    cpp/assign/resp.cpp
    cpp/assign/ring.cpp
    cpp/assign/tpacm4.cpp
    cpp/core/element_mass.cpp
    cpp/core/molecule.cpp
    cpp/core/template_ops.cpp
    cpp/forcefield/amber.cpp
    cpp/forcefield/amber_builtins.cpp
    cpp/forcefield/amber_names.cpp
    cpp/forcefield/amber_parser.cpp
    cpp/forcefield/amber_registry.cpp
    cpp/forcefield/amber_templates.cpp
    cpp/forcefield/charmm_loader.cpp
    cpp/forcefield/fep.cpp
    cpp/forcefield/gromacs_parser.cpp
    cpp/forcefield/nonamber.cpp
    cpp/forcefield/opls_loader.cpp
    cpp/forcefield/special_pairwise.cpp
    cpp/forcefield/special/gb.cpp
    cpp/io/coordinate.cpp
    cpp/io/mol2.cpp
    cpp/io/gro.cpp
    cpp/io/io_bundle.cpp
    cpp/io/pdb_hybrid36.cpp
    cpp/io/pdb.cpp
    cpp/io/pdb_reader.cpp
    cpp/io/pdb_writer.cpp
    cpp/io/psf.cpp
    cpp/io/sponge_coordinate.cpp
    cpp/io/sponge.cpp
    cpp/python/bindings.cpp
    cpp/python/bindings_assign.cpp
    cpp/python/bindings_core.cpp
    cpp/python/bindings_forcefield.cpp
    cpp/python/bindings_io.cpp
    cpp/solvation/solvation.cpp
    cpp/topology/topology.cpp
)

target_include_directories(_core PRIVATE cpp cpp/assign cpp/core cpp/forcefield cpp/io cpp/solvation cpp/topology)
if(TARGET HighFive::HighFive)
    target_link_libraries(_core PRIVATE HighFive::HighFive)
elseif(TARGET HighFive)
    target_link_libraries(_core PRIVATE HighFive)
else()
    message(FATAL_ERROR "HighFive target was not provided by find_package(HighFive)")
endif()
if(TARGET HDF5::HDF5)
    target_link_libraries(_core PRIVATE HDF5::HDF5)
else()
    target_include_directories(_core PRIVATE ${HDF5_INCLUDE_DIRS})
    target_link_libraries(_core PRIVATE ${HDF5_LIBRARIES})
endif()
target_compile_options(_core PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
)

option(XPONGECPP_SANITIZE "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
if (XPONGECPP_SANITIZE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(_core PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(_core PRIVATE -fsanitize=address,undefined)
endif()

install(TARGETS _core DESTINATION XpongeCPP)
