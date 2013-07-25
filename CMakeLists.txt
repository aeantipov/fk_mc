cmake_minimum_required(VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
project(ctint CXX)
set(CMAKE_BUILD_TYPE Release)

# find triqs
find_package(TRIQS REQUIRED)

# build fk mc library
add_library(fk_mc STATIC src/fk_mc.cpp)
include_directories(${TRIQS_INCLUDE_ALL})
target_link_libraries (fk_mc ${TRIQS_LIBRARY_ALL})
install(TARGETS fk_mc DESTINATION lib)

#add_subdirectory(test)
#enable_testing()
