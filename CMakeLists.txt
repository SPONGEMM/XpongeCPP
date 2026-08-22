cmake_minimum_required(VERSION 3.18)
project(XpongeCPP LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(PYBIND11_FINDPYTHON ON)

find_package(pybind11 CONFIG REQUIRED)
include(FetchContent)

# Build one private static HDF5 implementation on every platform.  In
# particular, this prevents the extension and h5py from loading two
# incompatible HDF5 shared libraries into the same Python process on Windows.
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries" FORCE)
set(BUILD_STATIC_LIBS ON CACHE BOOL "Build static libraries" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "Disable dependency tests" FORCE)
set(HDF5_BUILD_CPP_LIB OFF CACHE BOOL "Disable the HDF5 C++ API" FORCE)
set(HDF5_BUILD_EXAMPLES OFF CACHE BOOL "Disable HDF5 examples" FORCE)
set(HDF5_BUILD_FORTRAN OFF CACHE BOOL "Disable the HDF5 Fortran API" FORCE)
set(HDF5_BUILD_HL_LIB OFF CACHE BOOL "Disable the HDF5 high-level API" FORCE)
set(HDF5_BUILD_JAVA OFF CACHE BOOL "Disable the HDF5 Java API" FORCE)
set(HDF5_BUILD_TOOLS OFF CACHE BOOL "Disable HDF5 tools" FORCE)
set(HDF5_ENABLE_PARALLEL OFF CACHE BOOL "Disable parallel HDF5" FORCE)
set(HDF5_ENABLE_SZIP_SUPPORT OFF CACHE BOOL "Disable SZip support" FORCE)
set(HDF5_ENABLE_THREADSAFE OFF CACHE BOOL "Disable threaded HDF5" FORCE)
set(HDF5_ENABLE_Z_LIB_SUPPORT OFF CACHE BOOL "Disable zlib support" FORCE)
set(HDF5_NO_PACKAGES ON CACHE BOOL "Disable HDF5 packaging" FORCE)
FetchContent_Declare(
    hdf5
    URL https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_1.14.6.tar.gz
    URL_HASH SHA256=09ee1c671a87401a5201c06106650f62badeea5a3b3941e9b1e2e1e08317357f
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
)
FetchContent_MakeAvailable(hdf5)

set(HIGHFIVE_FIND_HDF5 OFF CACHE BOOL "Use XpongeCPP's private HDF5" FORCE)
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

pybind11_add_module(_core
    cpp/assign/amber_alternating.cpp
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
target_link_libraries(_core PRIVATE hdf5-static)
target_compile_options(_core PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
)
if(UNIX AND NOT APPLE)
    set_target_properties(_core PROPERTIES
        BUILD_RPATH_USE_ORIGIN TRUE
        INSTALL_RPATH "$ORIGIN"
        INSTALL_RPATH_USE_LINK_PATH FALSE
        INSTALL_REMOVE_ENVIRONMENT_RPATH TRUE
    )
endif()

option(XPONGECPP_SANITIZE "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
if (XPONGECPP_SANITIZE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(_core PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(_core PRIVATE -fsanitize=address,undefined)
endif()

install(TARGETS _core DESTINATION XpongeCPP)
install(
    DIRECTORY third_party/xponge_reference_forcefield/
    DESTINATION XpongeCPP/data/reference_forcefield
    PATTERN "__pycache__" EXCLUDE
    PATTERN "*.pyc" EXCLUDE
)
