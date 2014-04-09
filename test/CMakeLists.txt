include_directories(../solver .. ${TRIQS_INCLUDE_ALL})

set (tests_fk
#triqs_extra_test
lattice_test
config_test
binning_test
jackknife_test
#ietl_inverse01
ipr_test
mc_test01
)

set (other_tests
#eigen_test
hdf5
)

enable_testing()

foreach (test ${tests_fk})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    target_link_libraries(${test} fk_mc ${LINK_ALL})
    add_test(${test} ${test})
endforeach(test)

foreach (test ${other_tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    target_link_libraries(${test} ${LINK_ALL})
    add_test(${test} ${test})
endforeach(test)


