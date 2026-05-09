cmake_minimum_required(VERSION 3.18)
project(XpongeCPP LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(PYBIND11_FINDPYTHON ON)

find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_core
    cpp/core/molecule.cpp
    cpp/forcefield/amber.cpp
    cpp/io/mol2.cpp
    cpp/io/pdb.cpp
    cpp/io/sponge.cpp
    cpp/python/bindings.cpp
    cpp/solvation/solvation.cpp
    cpp/topology/topology.cpp
)

target_include_directories(_core PRIVATE cpp cpp/core cpp/forcefield cpp/io cpp/solvation cpp/topology)
target_compile_options(_core PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall -Wextra -Wpedantic>
)

option(XPONGECPP_SANITIZE "Build with AddressSanitizer and UndefinedBehaviorSanitizer" OFF)
if (XPONGECPP_SANITIZE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(_core PRIVATE -fsanitize=address,undefined -fno-omit-frame-pointer)
    target_link_options(_core PRIVATE -fsanitize=address,undefined)
endif()

install(TARGETS _core DESTINATION XpongeCPP)
