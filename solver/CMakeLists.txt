include_directories (.)

set (fk_mc_src
    lattice.hpp
    lattice/hypercubic.cpp
    lattice/triangular.cpp
    lattice/chain.cpp
    chebyshev.hpp
    configuration.hpp configuration.cpp
    moves.hpp moves.cpp
    moves_chebyshev.hpp moves_chebyshev.cpp
    measures/energy.cpp
    measures/spectrum.cpp
    measures/spectrum_history.cpp
    measures/focc_history.cpp
    measures/ipr.hpp
    measures/stiffness.hpp
    measures/eigenfunctions.cpp
    fk_mc.hpp fk_mc.cpp
)
add_library(fk_mc STATIC ${fk_mc_src})
target_link_libraries (fk_mc ${LINK_ALL})
install(TARGETS fk_mc DESTINATION lib)

install(DIRECTORY . DESTINATION include/fk_mc
        FILES_MATCHING PATTERN "*.hpp" PATTERN "*.hxx")
