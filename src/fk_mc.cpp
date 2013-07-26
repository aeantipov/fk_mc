#include "fk_mc.hpp"
#include "moves.hpp"
#include "measures.hpp"

#include <triqs/mc_tools/mc_generic.hpp>
#include <triqs/utility/callbacks.hpp>

namespace fk {

template <class lattice_t>
fk_mc<lattice_t>::fk_mc(lattice_t l):
    lattice(l)
{
}

template <class lattice_t>
void fk_mc<lattice_t>::solve(utility::parameters p)
{
    if (world.rank() == 0) std::cout << "Falicov-Kimball lattice Monte Carlo" << std::endl << std::endl;
    p.update(solve_defaults());

    mc_tools::mc_generic<double> mc(p);

    // Generate the configuration and cache the spectrum
    config_t config(lattice,p["U"],p["mu_c"],p["mu_f"]);
    config.randomize_f(mc.rng(),p["Nf_start"]);
    config.get_spectrum();

    mc.add_move(move_flip<lattice_t>(p["beta"], config, mc.rng()), "flip");
    mc.add_measure(dummy_measure(),"dummy");

      // run and collect results
    mc.start(1.0, triqs::utility::clock_callback(p["max_time"]));
    mc.collect_results(world);
}

template <class lattice_t>
 triqs::utility::parameter_defaults fk_mc<lattice_t>::solve_defaults() const {

  triqs::utility::parameter_defaults pdef;

  pdef.required
   ("beta", double(), "Inverse temperature")
   ("U", double(1.0), "FK U")
   ("N_Cycles", int(), "Number of QMC cycles")
   ;

  pdef.optional
   ("mu_c", double(0.5), "Chemical potential of c electrons")
   ("mu_f", double(0.5), "Chemical potential of f electrons")
   ("Nf_start", size_t(5), "Starting number of f-electrons")
   ("Length_Cycle", int(50), "Length of a single QMC cycle")
   ("N_Warmup_Cycles", int(5000), "Number of cycles for thermalization")
   ("Random_Seed", int(34788+928374*world.rank()), "Seed for random number generator")
   ("Random_Generator_Name", std::string(""), "Name of random number generator")
   ("max_time",int(600), "Maximum running time")
   ;

  return pdef;
 }

/* Explicit instantiation */
template class fk_mc<square_lattice_traits<2>>; // Square lattice, 2d

} // end of namespace FK
