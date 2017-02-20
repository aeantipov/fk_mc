# Build executables

#find_package (TCLAP)
#message(STATUS "TCLAP includes: " ${TCLAP_INCLUDE_DIRS} )

link_libraries(fk_mc ${LINK_ALL})
include_directories(../solver ${TRIQS_INCLUDE_ALL} ${TCLAP_INCLUDE_DIRS})

#set (progs 
#fk_mc_exec
#)

#lattices
foreach (lattice ${lattices})
#    message(STATUS "Lattice: ${lattice}")
#    set (prog_src data_save.hpp data_save.hxx data_load.hpp fk_mc_exec.cpp)
#    set (prog fk_mc_${lattice})
#    add_executable(${prog} ${prog_src})
#    set_target_properties(${prog} PROPERTIES COMPILE_DEFINITIONS "LATTICE_${lattice}")
#    install ( TARGETS ${prog} DESTINATION bin )

    set (prog_src data_save.hpp data_save.hxx data_load.hpp fk_mc_kwant.cpp)
    set (prog fk_kwant_${lattice})
    add_executable(${prog} ${prog_src})
    set_target_properties(${prog} PROPERTIES COMPILE_DEFINITIONS "LATTICE_${lattice}") 
    install ( TARGETS ${prog} DESTINATION bin )
endforeach(lattice)



