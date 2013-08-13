cmake_minimum_required(VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
project(ctint CXX)
set(CMAKE_BUILD_TYPE Release)

# find triqs
find_package(TRIQS REQUIRED)

find_package (Arpack)
message(STATUS "Arpack libs: " ${ARPACK_LIBRARIES})
find_package (Eigen3 3.1)
message(STATUS "Eigen3 includes: " ${EIGEN3_INCLUDE_DIR} )

set (LINK_ALL ${TRIQS_LIBRARY_ALL} ${ARPACK_LIBRARIES})
include_directories(${TRIQS_INCLUDE_ALL} ${EIGEN3_INCLUDE_DIR})

# build fk mc library
add_library(fk_mc STATIC solver/lattice_traits.cpp solver/fk_mc.cpp)
target_link_libraries (fk_mc ${LINK_ALL})
install(TARGETS fk_mc DESTINATION lib)

# build statistics
#add_subdirectory(statistics)

# build executable
option(Execs "Enable executables" ON)
if (Execs)
    add_subdirectory(prog)
endif (Execs)

#build tests
option(Testing "Enable testing" OFF)
if (Testing)
    message(STATUS "Enabled test building and running")
    add_subdirectory(test)
    enable_testing()
endif (Testing)
