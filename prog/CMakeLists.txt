# Build executables

find_package (TCLAP)
message(STATUS "TCLAP includes: " ${TCLAP_INCLUDE_DIRS} )

link_libraries(fk_mc ${TRIQS_LIBRARY_ALL})
include_directories(../solver ../statistics ${TRIQS_INCLUDE_ALL} ${TCLAP_INCLUDE_DIRS})

set (progs 
fk_mc_exec
)

foreach (prog ${progs})
    set(prog_src ${prog}.cpp)
    add_executable(${prog} ${prog_src})
    install ( TARGETS ${prog} DESTINATION bin )
endforeach(prog)


