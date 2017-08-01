# Build executables
foreach (lattice ${lattices})
    message(STATUS "Lattice: ${lattice}")
    set (prog_src data_save.hpp data_save.hxx data_load.hpp fk_mc_exec.cpp)
    set (prog fk_mc_${lattice})
    add_executable(${prog} ${prog_src})
    set_target_properties(${prog} PROPERTIES COMPILE_DEFINITIONS "LATTICE_${lattice}") 
    target_include_directories(${prog} PUBLIC $<TARGET_PROPERTY:${PROJECT_NAME},INTERFACE_INCLUDE_DIRECTORIES>)
    target_link_libraries(${prog} PUBLIC ${PROJECT_NAME} ${${PROJECT_NAME}_DEPENDS})
    install ( TARGETS ${prog} DESTINATION bin )
endforeach(lattice)



