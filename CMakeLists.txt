cmake_minimum_required(VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
project(ctint CXX)
set(CMAKE_BUILD_TYPE Release)

# find triqs
find_package(TRIQS REQUIRED)

include_directories(${TRIQS_INCLUDE_ALL})

# build fk mc library
add_library(fk_mc STATIC src/lattice_traits.cpp src/fk_mc.cpp)
target_link_libraries (fk_mc ${TRIQS_LIBRARY_ALL})
install(TARGETS fk_mc DESTINATION lib)

#build tests
enable_testing()
link_libraries(fk_mc ${TRIQS_LIBRARY_ALL})
set (tests 
lattice_traits_test
config_test
mc_test01
)
include_directories(./src)
foreach (test ${tests})
    set(test_src test/${test}.cpp)
    add_executable(test/${test} ${test_src})
    
    add_test(${test} test/${test})
endforeach(test)

