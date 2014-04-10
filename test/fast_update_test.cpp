#include <boost/mpi/environment.hpp>
#include <chrono>
#include <Eigen/Sparse>

#include "fk_mc.hpp"
#include "moves.hpp"

using namespace fk;

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
  try {

    size_t L = 11;

    double U = 1.0;
    double mu = U/2;
    double e_f = 0.0;
    double t = -1.0;
    double T = 0.15;
    double beta = 1.0/T;

    triqs::mc_tools::random_generator r("mt19937", 32167);

    typedef hypercubic_lattice<2> lattice_t;
    lattice_t lattice(L);
    lattice.fill(-1.0);

    configuration_t config(lattice, 1.0, U, mu, mu+e_f);
    config.randomize_f(r,L*L/2);
    config.calc_hamiltonian();
    config.calc_full_spectrum(true);

    configuration_t config1(config);
    config1 = config;

    move_addremove move1(beta, config, r);
    move1.attempt();
    move1.accept();

    MY_DEBUG("Old config: " << config1.f_config.transpose());
    MY_DEBUG("New config: " << config.f_config.transpose());
    MY_DEBUG("Old spectrum: " << config1.cached_spectrum.transpose());
    MY_DEBUG("New spectrum: " << config.cached_spectrum.transpose() + 0.0*(config1.cached_spectrum[0] - config.cached_spectrum[0]));
    
    typedef configuration_t::dense_m dense_m;
    if (0==1) { 
        dense_m a = config.hamilt;
        dense_m diag(a); a.setZero(); a.diagonal().setOnes();
        MY_DEBUG((a + diag*2.) * (a - diag*2.).inverse())
        
        }
   }
  catch(triqs::runtime_error const & e) { std::cout  << "exception "<< e.what() << std::endl;}
  return 0;

    return EXIT_SUCCESS;
}
