# Build executables

link_libraries(fk_mc ${LINK_ALL})

#lattices
foreach (lattice ${lattices})
    message(STATUS "Lattice: ${lattice}")
    set (prog_src data_save.hpp data_save.hxx data_load.hpp fk_mc_exec.cpp)
    set (prog fk_mc_${lattice})
    add_executable(${prog} ${prog_src})
    set_target_properties(${prog} PROPERTIES COMPILE_DEFINITIONS "LATTICE_${lattice}")
    install ( TARGETS ${prog} DESTINATION bin )
endforeach(lattice)



