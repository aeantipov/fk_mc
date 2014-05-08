#include "fk_mc.hpp"
#include "measures/fsusc_fft.hpp"
#include "measures/fsusc0pi.hpp"
#include <boost/mpi/environment.hpp>
#include <chrono>


using namespace fk;

//complex_type glocal_imfreq(ComplexType z, std::vector<
//calc_eigenvectors_);
struct dummy_measure_a
{
    dummy_measure_a(){};
    void accumulate(double sign){std::cout << "dummy_measure_a" << std::endl;};
    void collect_results(boost::mpi::communicator const &c){};
};

struct dummy_measure_b
{
    dummy_measure_b(){};
    void accumulate(double sign){std::cout << "dummy_measure_b" << std::endl;};
    void collect_results(boost::mpi::communicator const &c){};
};

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
  try {
    size_t L = 4;

    double U = 1.0;
    double t = -1.0;
    double T = 0.1;
    double beta = 1.0/T;

    triangular_lattice lattice(L);
    lattice.fill(-1.0,0);


    triqs::utility::parameters p;
    p["U"] = U;
    p["mu_c"] = U/2; p["mu_f"] = U/2;
    p["beta"] = beta;
    p["Nf_start"] = L*L/2;
    p["random_seed"] = 34788;
    p["verbosity"] = 3;
    p["length_cycle"] = 1; 
    p["n_warmup_cycles"] = 1;
    p["n_cycles"] = 40;
    p["max_time"]=5;
    p["measure_history"] = true;

    FKDEBUG(p);
    fk_mc mc(lattice,p);

    std::vector<double> n0, npi, n0_n0, npi_npi;

    //mc.add_measure(measure_fsusc<triangular_lattice>(beta, lattice, mc.config, {lattice.get_bzpoint({0, 0}), lattice.get_bzpoint({PI, PI})}, mc.observables.nq_history, mc.observables.fsuscq_history), "fsusc");
    mc.add_measure(measure_nf0pi<triangular_lattice>(mc.config, lattice, mc.observables.nf0, mc.observables.nfpi), "nf0pi");
    //mc.add_measure(dummy_measure,"a");
    mc.add_measure(dummy_measure_a(),"b");
    mc.add_measure(dummy_measure_b(),"a");
    mc.solve();

    const auto &spectrum_history = mc.observables.spectrum_history;
    for (auto x: spectrum_history[4]) std::cout << x << " " << std::flush; std::cout << std::endl;

   /* const auto &nq_history = mc.observables.fsuscq_history;
    const auto &fsuscq_history = mc.observables.fsuscq_history;
    for (auto x : fsuscq_history[0]) std::cout << x << " " << std::flush; std::cout << std::endl;
    for (auto x : fsuscq_history[1]) std::cout << x << " " << std::flush; std::cout << std::endl;
    */

    //auto e = mc

  }
  catch(triqs::runtime_error const & e) { std::cout  << "exception "<< e.what() << std::endl;}
  catch(std::exception const & e) { std::cout  << "exception "<< e.what() << std::endl;}
  return 0;

    return EXIT_SUCCESS;
}
