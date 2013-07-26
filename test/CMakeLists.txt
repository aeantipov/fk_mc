link_libraries(fk_mc ${TRIQS_LIBRARY_ALL})
include_directories(../src ${TRIQS_INCLUDE_ALL})

set (tests 
lattice_traits_test
config_test
mc_test01
)

enable_testing()

foreach (test ${tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    add_test(${test} ${test})
endforeach(test)


