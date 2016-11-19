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
#include "measures/stiffness.hpp"
#include "measures/eigenfunctions.hpp"

//#include <triqs/utility/callbacks.hpp>

namespace fk {

//template <typename L>
//parameters_t _update_def(triqs::utility::parameters p){p.update(fk_mc<L>::solve_defaults()); return p;}

template <typename L>
//fk_mc<L>::fk_mc(lattice_type l, triqs::utility::parameters p1, bool randomize_config):
fk_mc<L>::fk_mc(parameters_t const &p, int rank):
    base(p, rank), 
    p(p)
{
}

template <typename L>
//fk_mc<L>::fk_mc(lattice_type l, triqs::utility::parameters p1, bool randomize_config):
void fk_mc<L>::initialize(lattice_type l, bool randomize_config, std::vector<double> wgrid_conductivity)
    //mc(p) 
{
    lattice_ptr = std::make_shared<lattice_type>(l);
    config_ptr = std::make_shared<configuration_t> (configuration_t(*lattice_ptr,p["beta"],p["U"],p["mu_c"],p["mu_f"]));
    std::cout << "\tRandom seed for proc " << comm.rank() << " : " << p["seed"] << std::endl;
    if (randomize_config) config_ptr->randomize_f(this->rng(),p["Nf_start"]);

    configuration_t& config = *config_ptr;
    lattice_type& lattice = *lattice_ptr;

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

    size_t max_bins = p["nsweeps"];
    observables.reserve(max_bins);

    bool calc_spectrum = !cheb_move; 
    this->add_measure(measure_nf0pi<lattice_type>(config, lattice, observables.nf0, observables.nfpi), "nf0pi");
    if (p["measure_history"]) { if (!comm.rank()) std::cout << "Saving history" << std::endl; };
    if (p["measure_history"] && p["measure_ipr"]) {
        if (!comm.rank()) std::cout << "Measuring ipr" << std::endl;
        this->add_measure(measure_ipr<lattice_type>(config, lattice, observables.ipr_history),"ipr");
        }
    if (p["measure_stiffness"]) {
        if (!comm.rank()) std::cout << "Measuring stiffness" << std::endl;
        this->add_measure(measure_stiffness<lattice_type>(config, lattice, observables.stiffness, observables.cond_history, wgrid_conductivity, p["cond_offset"]),"stiffness");
        };

    if (p["measure_eigenfunctions"]) {
        if (!comm.rank()) std::cout << "Measuring eigenfunctions" << std::endl;
        this->add_measure(measure_eigenfunctions(config,observables.eigenfunctions_history), "eigenfunctions");
        };

    calc_spectrum = calc_spectrum || p["measure_ipr"];

    if (!cheb_move || calc_spectrum) {
        this->add_measure(measure_energy(beta,config,observables.energies, observables.d2energies, observables.c_energies), "energy");
        this->add_measure(measure_spectrum(config,observables.spectrum), "spectrum");
        if (p["measure_history"]) 
            this->add_measure(measure_spectrum_history(config,observables.spectrum_history), "spectrum_history");
        };
    if (p["measure_history"]) {
        this->add_measure(measure_focc(config,observables.focc_history), "focc_history");
        };
}

/*
template <typename L>
void fk_mc<L>::solve(std::vector<double> wgrid_conductivity)
{
    if (int(p["n_cycles"]) == 0) return;
    if (comm.rank() == 0) std::cout << "Running MC..." << std::endl << std::endl;

          // run and collect results
    this->start(1.0, [](){return false;}); 
    //this->start(1.0, triqs::utility::clock_callback(p["max_time"]));
    this->collect_results(comm);
}
*/

/*template <typename L>
 triqs::utility::parameter_defaults fk_mc<L>::solve_defaults() {

  triqs::utility::parameter_defaults pdef;

  pdef.required
   ("beta", double(), "Inverse temperature")
   .required("U", double(1.0), "FK U")
   .required("n_cycles", int(), "Number of QMC cycles")
   .required("mu_c", double(0.5), "Chemical potential of c electrons")
   .required("mu_f", double(0.5), "Chemical potential of f electrons")
   ;

  pdef.optional
   ("mc_flip", double(0.0), "Make flip moves")
   .optional("mc_add_remove", double(1.0), "Make add/remove moves")
   .optional("mc_reshuffle", double(0.0), "Make reshuffle moves")
   .optional("cheb_moves", bool(false), "Allow moves using Chebyshev sampling")
   .optional("cheb_prefactor", double(2.2), "Prefactor for number of Chebyshev polynomials = #ln(Volume)")
   .optional("measure_history", bool(true), "Measure the history")
   .optional("random_name", std::string(""), "Name of random number generator")
   .optional("Nf_start", size_t(5), "Starting number of f-electrons")
   .optional("length_cycle", int(50), "Length of a single QMC cycle")
   .optional("n_warmup_cycles", int(5000), "Number of cycles for thermalization")
   .optional("random_seed", int(34788), "Seed for random number generator")
   .optional("max_time",int(600000), "Maximum running time")
   .optional("measure_ipr", bool(false), "Measure inverse participation ratio")
   .optional("measure_eigenfunctions", bool(false), "Measure eigenfunctions")
   .optional("cond_offset", double(0.05), "dos offset from the real axis")
   ;

  return pdef;
 }
*/

template <typename L>
parameters_t &fk_mc<L>::define_parameters(parameters_type &p) {
   base::define_parameters(p);

   p.define<double> ("beta", double(10.0), "Inverse temperature")
   .define<double> ("U", double(1.0), "FK U")
   //.define<int>("n_cycles", int(), "Number of QMC cycles")
   .define<double>("mu_c", double(0.5), "Chemical potential of c electrons")
   .define<double>("mu_f", double(0.5), "Chemical potential of f electrons")
   ;

   p.define<double>("mc_flip", double(0.0), "Make flip moves")
   .define<double>("mc_add_remove", double(1.0), "Make add/remove moves")
   .define<double>("mc_reshuffle", double(0.0), "Make reshuffle moves")
   .define<bool>("cheb_moves", bool(false), "Allow moves using Chebyshev sampling")
   .define<double>("cheb_prefactor", double(2.2), "Prefactor for number of Chebyshev polynomials = #ln(Volume)")
   .define<bool>("measure_history", bool(true), "Measure the history")
   //.optional("random_name", std::string(""), "Name of random number generator")
   .define<int>("Nf_start", size_t(5), "Starting number of f-electrons")
   //.define<int>("length_cycle", int(50), "Length of a single QMC cycle")
   //.define<int>("n_warmup_cycles", int(5000), "Number of cycles for thermalization")
   .define<int>("seed", int(std::random_device()()), "Seed for random number generator")
   .define<bool>("measure_ipr", bool(false), "Measure inverse participation ratio")
   .define<bool>("measure_eigenfunctions", bool(false), "Measure eigenfunctions")
   .define<double>("cond_offset", double(0.05), "dos offset from the real axis")
   .define<bool>("measure_stiffness", bool(false), "Measure stiffness/conductivity")
   ;

  p["SEED"] = p["seed"];

  return p;
 }



} // end of namespace FK
