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

add_definitions("-Wno-deprecated-register -Wno-mismatched-new-delete")

include_directories(${CMAKE_SOURCE_DIR}/include/fk_mc ${CMAKE_BINARY_DIR}/include/fk_mc)
add_this_package(${fk_mc_src})

# Dependencies
find_package (Arpack)
message(STATUS "Arpack libraries: " ${ARPACK_LIBRARIES} )
link_libraries(${ARPACK_LIBRARIES})

add_alpscore()
add_eigen3()
add_gftools()
add_boost(mpi)
add_fftw3()

#configure compile definitions file
configure_file(${CMAKE_SOURCE_DIR}/include/fk_mc/definitions.hpp.in ${CMAKE_BINARY_DIR}/include/fk_mc/definitions.hpp)

# build fk mc library
#add_subdirectory(solver)

# build executable
option(Execs "Enable executables" ON)
if (Execs)
    list(APPEND all_lattices cubic1d cubic2d cubic3d)
    list(LENGTH LATTICES lattice_len)
    if (${lattice_len} EQUAL 0) 
        set(LATTICES ${all_lattices})
    endif()
    set (lattices ${LATTICES} CACHE STRING "Lattices for FK MC calculation" FORCE)
    add_subdirectory(prog)
endif (Execs)

option(Benchmark "Enable benchmark" OFF)
option(Testing "Enable testing" OFF)

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
