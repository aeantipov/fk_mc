link_libraries(fk_mc ${TRIQS_LIBRARY_ALL})
include_directories(../solver .. ${TRIQS_INCLUDE_ALL})

set (tests 
triqs_extra_test
lattice_traits_test
config_test
binning_test
jackknife_test
#ietl_inverse01
#mc_test01
)

enable_testing()

foreach (test ${tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    add_test(${test} ${test})
endforeach(test)


