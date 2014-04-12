cmake_minimum_required(VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
project(fk_mc CXX)

# Toggle debug mode on/off
option (Debug "Debug mode" OFF)
if (Debug)
    message (STATUS "Debug mode")
    set(CMAKE_BUILD_TYPE Debug)
    set(FK_MC_DEBUG true)
    set(Testing true FORCE)
else (Debug)
    message (STATUS "Release mode")
    set(CMAKE_BUILD_TYPE Release)
endif(Debug)

# RPATH settings (see http://www.cmake.org/Wiki/CMake_RPATH_handling) 
SET(CMAKE_SKIP_BUILD_RPATH  FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE) 
SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
LIST(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES
    "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
IF("${isSystemDir}" STREQUAL "-1")
    SET(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
ENDIF("${isSystemDir}" STREQUAL "-1")

# find triqs
find_package(TRIQS REQUIRED)
message(STATUS ${TRIQS_CXX_DEFINITIONS})
add_definitions ( ${TRIQS_CXX_DEFINITIONS})

# find Eigen 
find_package (Eigen3 3.1)
message(STATUS "Eigen3 includes: " ${EIGEN3_INCLUDE_DIR} )
find_package (Arpack)
message(STATUS "Arpack libraries: " ${ARPACK_LIBRARIES} )
#find FFTW
#find_package (FFTW)
#message(STATUS "FFTW includes: " ${FFTW_INCLUDE_DIR} )

set (LINK_ALL ${TRIQS_LIBRARY_ALL} ${ARPACK_LIBRARIES} ${FFTW_LIBRARIES})
include_directories(${TRIQS_INCLUDE_ALL} ${EIGEN3_INCLUDE_DIR} ${FFTW_INCLUDE_DIR} ${CMAKE_BINARY_DIR}/solver)


#configure compile definitions file
configure_file(${CMAKE_SOURCE_DIR}/solver/definitions.hpp.in ${CMAKE_BINARY_DIR}/solver/definitions.hpp)

# build fk mc library
add_subdirectory(solver)

# build executable
option(Execs "Enable executables" ON)
if (Execs)
    set (lattices ${LATTICES} CACHE STRING "Lattices for FK MC calculation" FORCE)
    add_subdirectory(prog)
endif (Execs)

#build tests
option(Testing "Enable testing" OFF)
if (Testing)
    message(STATUS "Enabled test building and running")
    add_subdirectory(test)
    enable_testing()
endif (Testing)
