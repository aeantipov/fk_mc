include_directories(../solver .. ${TRIQS_INCLUDE_ALL})

set (tests_fk
lattice_test
honeycomb_test
#chain_test
config_test
#fast_update_test
binning_test
jackknife_test
ipr_test
stiffness_test
graph_lattice_test
#polarized_test
#mc_test01
#saveload_test
)

set (other_tests
#eigen_test
#hdf5
)

enable_testing()

foreach (test ${tests_fk})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    target_link_libraries(${test} fk_mc gtest_main ${LINK_ALL})
    add_test(${test} ${test})
endforeach(test)

foreach (test ${other_tests})
    set(test_src ${test}.cpp)
    add_executable(${test} ${test_src})
    target_link_libraries(${test} gtest_main ${LINK_ALL})
    add_test(${test} ${test})
endforeach(test)


