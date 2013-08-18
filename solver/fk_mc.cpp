#include "fk_mc.hpp"
#include "moves.hpp"
#include "measure_energy.hpp"
#include "measure_d2energy.hpp"

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
    if (world.rank() == 0) std::cout << "Running MC..." << std::endl << std::endl;
    p.update(solve_defaults());
    INFO("\tRandom seed for proc " << world.rank() << " : " << p["Random_Seed"]);

    mc_tools::mc_generic<double> mc(p);

    // Generate the configuration and cache the spectrum
    config_t config(lattice,p["U"],p["mu_c"],p["mu_f"]);
    config.randomize_f(mc.rng(),p["Nf_start"]);
    config.calc_hamiltonian();
    config.calc_spectrum();

    double beta = p["beta"];
    if (double(p["mc_flip"])>std::numeric_limits<double>::epsilon()) 
        mc.add_move(move_flip<config_t>(beta, config, mc.rng()), "flip", p["mc_flip"]);
    if (double(p["mc_add_remove"])>std::numeric_limits<double>::epsilon()) 
        mc.add_move(move_addremove<config_t>(beta, config, mc.rng()), "add_remove", p["mc_add_remove"]);
    if (double(p["mc_reshuffle"])>std::numeric_limits<double>::epsilon()) 
        mc.add_move(move_randomize<config_t>(beta, config, mc.rng()), "reshuffle", p["mc_reshuffle"]);

    size_t max_bins = p["N_Cycles"];
    observables.energies.reserve(max_bins);
    mc.add_measure(measure_energy<config_t>(beta,config,observables.energies), "energy");
    mc.add_measure(measure_d2energy<config_t>(beta,config,observables.d2energies), "d2energy");

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
   ("mu_c", double(0.5), "Chemical potential of c electrons")
   ("mu_f", double(0.5), "Chemical potential of f electrons")
   ;

  pdef.optional
   ("mc_flip", double(1.0), "Make flip moves")
   ("mc_add_remove", double(1.0), "Make add/remove moves")
   ("mc_reshuffle", double(1.0), "Make reshuffle moves")
   ("Nf_start", size_t(5), "Starting number of f-electrons")
   ("Length_Cycle", int(50), "Length of a single QMC cycle")
   ("N_Warmup_Cycles", int(5000), "Number of cycles for thermalization")
   ("Random_Seed", int(34788), "Seed for random number generator")
   ("Random_Generator_Name", std::string(""), "Name of random number generator")
   ("max_time",int(600), "Maximum running time")
   ;

  return pdef;
 }

/* Explicit instantiation */
//template class fk_mc<square_lattice_traits<2>>; // Square lattice, 2d
template class fk_mc<triangular_lattice_traits>; // Triangular lattice, 2d

} // end of namespace FK
