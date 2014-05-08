include_directories(../solver .. ${TRIQS_INCLUDE_ALL})

set (benchmarks
fast_update
)

foreach (benchmark ${benchmarks})
    set(benchmark_src ${benchmark}.cpp)
    add_executable(${benchmark} ${benchmark_src})
    set_target_properties(${benchmark} PROPERTIES COMPILE_FLAGS "-DNDEBUG -O3")
    target_link_libraries(${benchmark} fk_mc gtest_main ${LINK_ALL})
endforeach(benchmark)


