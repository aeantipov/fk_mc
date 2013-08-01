cmake_minimum_required(VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
project(ctint CXX)
set(CMAKE_BUILD_TYPE Release)

# find triqs
find_package(TRIQS REQUIRED)

include_directories(${TRIQS_INCLUDE_ALL})

# build fk mc library
add_library(fk_mc STATIC solver/lattice_traits.cpp solver/fk_mc.cpp)
target_link_libraries (fk_mc ${TRIQS_LIBRARY_ALL})
install(TARGETS fk_mc DESTINATION lib)

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
