include_directories (.)

set (fk_mc_src
    lattice.cpp
    configuration.cpp
    moves.cpp
    measures/energy.cpp
    measures/spectrum.cpp
    measures/spectrum_history.cpp
    measures/focc_history.cpp
    fk_mc.cpp
)
add_library(fk_mc STATIC ${fk_mc_src})
target_link_libraries (fk_mc ${LINK_ALL})
install(TARGETS fk_mc DESTINATION lib)

#install(DIRECTORY . DESTINATION include/fk_mc
#          FILES_MATCHING PATTERN "*.h")
