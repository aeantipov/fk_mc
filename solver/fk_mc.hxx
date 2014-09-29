#pragma once

#include "fk_mc.hpp"
#include "moves.hpp"
#include "moves_chebyshev.hpp"
#include "measures/energy.hpp"
#include "measures/spectrum.hpp"
#include "measures/spectrum_history.hpp"
#include "measures/focc_history.hpp"
#include "measures/fsusc0pi.hpp"
#include "measures/ipr.hpp"

#include <triqs/utility/callbacks.hpp>

namespace fk {

template <typename L>
triqs::utility::parameters _update_def(triqs::utility::parameters p){p.update(fk_mc<L>::solve_defaults()); return p;}

template <typename L>
fk_mc<L>::fk_mc(lattice_type l, utility::parameters p1, bool randomize_config):
    base(_update_def<L>(p1)),
    p(_update_def<L>(p1)),
    lattice(l),
    //cheb_eval(chebyshev_eval(p["cheb_size"], p["cheb_grid_size"])),
    config(lattice,p["beta"],p["U"],p["mu_c"],p["mu_f"])
    //mc(p) 
{
    INFO("\tRandom seed for proc " << comm.rank() << " : " << p["random_seed"]);
    if (randomize_config) config.randomize_f(this->rng(),p["Nf_start"]);
}

template <typename L>
void fk_mc<L>::solve()
{
    if (int(p["n_cycles"]) == 0) return;
    if (comm.rank() == 0) std::cout << "Running MC..." << std::endl << std::endl;

    // Generate the configuration_t and cache the spectrum
    double beta = p["beta"];
    config.calc_hamiltonian();

    std::unique_ptr<chebyshev::chebyshev_eval> cheb_ptr;

    bool cheb_move = p["cheb_moves"];
    if (cheb_move) {
        int cheb_size = int(std::log(lattice.get_msize()) * double(p["cheb_prefactor"]));
        cheb_size+=cheb_size%2;
        size_t ngrid_points = std::max(cheb_size*2,10);
        cheb_ptr.reset(new chebyshev::chebyshev_eval(cheb_size, ngrid_points));
    }
        

    if (double(p["mc_flip"])>std::numeric_limits<double>::epsilon()) { 
        if (!cheb_move) this->add_move(move_flip(beta, config, this->rng()), "flip", p["mc_flip"]); 
                   else this->add_move(chebyshev::move_flip(beta, config, *cheb_ptr, this->rng()), "flip", p["mc_flip"]); 
        };
    if (double(p["mc_add_remove"])>std::numeric_limits<double>::epsilon()) { 
        if (!cheb_move) this->add_move(move_addremove(beta, config, this->rng()), "add_remove", p["mc_add_remove"]);
                   else this->add_move(chebyshev::move_addremove(beta, config, *cheb_ptr, this->rng()), "add_remove", p["mc_add_remove"]);
        };
    if (double(p["mc_reshuffle"])>std::numeric_limits<double>::epsilon()) { 
        if (!cheb_move) this->add_move(move_randomize(beta, config, this->rng()),  "reshuffle", p["mc_reshuffle"]);
                   else this->add_move(chebyshev::move_randomize(beta, config, *cheb_ptr, this->rng()), "reshuffle", p["mc_reshuffle"]);
        };

    size_t max_bins = p["n_cycles"];
    observables.reserve(max_bins);

    bool calc_spectrum = !cheb_move; 
    this->add_measure(measure_nf0pi<lattice_type>(config, lattice, observables.nf0, observables.nfpi), "nf0pi");
    if (p["measure_history"] && p["measure_ipr"]) {
        if (!comm.rank()) std::cout << "Measuring ipr" << std::endl;
        this->add_measure(measure_ipr<lattice_type>(config, lattice, observables.ipr_history),"ipr");
        }

    calc_spectrum = calc_spectrum || p["measure_ipr"];

    if (!cheb_move || calc_spectrum) {
        this->add_measure(measure_energy(beta,config,observables.energies, observables.d2energies), "energy");
        this->add_measure(measure_spectrum(config,observables.spectrum), "spectrum");
        if (p["measure_history"]) 
            this->add_measure(measure_spectrum_history(config,observables.spectrum_history), "spectrum_history");
        };
    if (p["measure_history"]) {
        this->add_measure(measure_focc(config,observables.focc_history), "focc_history");
        };

      // run and collect results
    this->start(1.0, [](){return false;}); 
    //this->start(1.0, triqs::utility::clock_callback(p["max_time"]));
    this->collect_results(comm);
}

template <typename L>
 triqs::utility::parameter_defaults fk_mc<L>::solve_defaults() {

  triqs::utility::parameter_defaults pdef;

  pdef.required
   ("beta", double(), "Inverse temperature")
   ("U", double(1.0), "FK U")
   ("n_cycles", int(), "Number of QMC cycles")
   ("mu_c", double(0.5), "Chemical potential of c electrons")
   ("mu_f", double(0.5), "Chemical potential of f electrons")
   ;

  pdef.optional
   ("mc_flip", double(0.0), "Make flip moves")
   ("mc_add_remove", double(1.0), "Make add/remove moves")
   ("mc_reshuffle", double(0.0), "Make reshuffle moves")
   ("cheb_moves", bool(false), "Allow moves using Chebyshev sampling")
   ("cheb_prefactor", double(2.2), "Prefactor for number of Chebyshev polynomials = #ln(Volume)")
   ("measure_history", bool(true), "Measure the history")
   ("random_name", std::string(""), "Name of random number generator")
   ("Nf_start", size_t(5), "Starting number of f-electrons")
   ("length_cycle", int(50), "Length of a single QMC cycle")
   ("n_warmup_cycles", int(5000), "Number of cycles for thermalization")
   ("random_seed", int(34788), "Seed for random number generator")
   ("max_time",int(600000), "Maximum running time")
   ("measure_ipr", bool(false), "Measure inverse participation ratio")
   ;

  return pdef;
 }

} // end of namespace FK
