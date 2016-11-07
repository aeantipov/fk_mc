cmake_minimum_required(VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
project(fk_mc CXX)

set(FK_MC_VERSION 0.1)

# check ../cmake/BoldHybCommondModuleDefinitions.cmake for all of the macro/functions below
include(CommonDefs)

# Disable build in source
no_source_builds()
# Rpath to absoulte path
fix_rpath()
# Workarounds for bugs in compilers
compiler_workarounds()
# C++11
option(AutoSetCXX11 "Auto set c++11 flags" ON)
if (AutoSetCXX11)
    set_cxx11()
endif (AutoSetCXX11)
# Print build type
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set (fk_mc_src
    mc_metropolis
    lattice/hypercubic
    lattice/triangular
    lattice/chain
    lattice/honeycomb
    configuration
    moves
    moves_chebyshev
    measures/energy
    measures/spectrum
    measures/spectrum_history
    measures/focc_history
    measures/eigenfunctions
    fk_mc
)

include_directories(${CMAKE_SOURCE_DIR}/include/fk_mc ${CMAKE_BINARY_DIR}/include/fk_mc)
add_this_package(${fk_mc_src})

# Dependencies
add_alpscore()
add_eigen3()
add_gftools()


# find triqs
#find_package(TRIQS REQUIRED)
#message(STATUS ${TRIQS_CXX_DEFINITIONS})
#add_definitions ( ${TRIQS_CXX_DEFINITIONS})

# find Eigen 
#find_package (Eigen3 3.1)
#message(STATUS "Eigen3 includes: " ${EIGEN3_INCLUDE_DIR} )
find_package (Arpack)
message(STATUS "Arpack libraries: " ${ARPACK_LIBRARIES} )
link_libraries(${ARPACK_LIBRARIES})
#find FFTW
#find_package (FFTW)
#message(STATUS "FFTW includes: " ${FFTW_INCLUDE_DIR} )

#set (LINK_ALL ${TRIQS_LIBRARY_ALL} ${ARPACK_LIBRARIES} ${FFTW_LIBRARIES})
#include_directories(${TRIQS_INCLUDE_ALL} ${EIGEN3_INCLUDE_DIR} ${FFTW_INCLUDE_DIR} ${CMAKE_BINARY_DIR}/solver)
#target_include_directories(${PROJEC

add_boost(mpi)

#configure compile definitions file
configure_file(${CMAKE_SOURCE_DIR}/include/fk_mc/definitions.hpp.in ${CMAKE_BINARY_DIR}/include/fk_mc/definitions.hpp)

# build fk mc library
#add_subdirectory(solver)

# build executable
option(Execs "Enable executables" ON)
if (Execs)
    list(APPEND all_lattices cubic2d triangular cubic3d honeycomb)
    list(LENGTH LATTICES lattice_len)
    if (${lattice_len} EQUAL 0) 
        set(LATTICES ${all_lattices})
    endif()
    set (lattices ${LATTICES} CACHE STRING "Lattices for FK MC calculation" FORCE)
    add_subdirectory(prog)
endif (Execs)

option(Benchmark "Enable benchmark" ON)
option(Testing "Enable testing" ON)

if (Benchmark OR Testing)
    include(EnableGtests) #defined in common/cmake
    #include(UseGtest)
    #add_subdirectory(${gtest_ROOT} ${PROJECT_BINARY_DIR}/gtest)
    include_directories(${gtest_INCLUDE_DIR})
endif()

if (Benchmark)
    message(STATUS "Building benchmark")
    add_subdirectory(benchmark)
endif (Benchmark)

#build tests
if (Testing)
    enable_testing()
    add_subdirectory(test)
endif(Testing)


### Additional installation files
configure_file("${CMAKE_SOURCE_DIR}/fk_mc.lmod.in" "${CMAKE_BINARY_DIR}/${FK_MC_VERSION}")
