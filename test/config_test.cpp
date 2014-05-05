#include "fk_mc.hpp"
#include <boost/mpi/environment.hpp>
#include <chrono>

using namespace fk;

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
  try {

    size_t L = 4;

    double U = 1.0;
    double mu = U/2;
    double e_f = 0.0;
    double t = -1.0;
    double T = 0.15;
    double beta = 1.0/T;

    hypercubic_lattice<2> lattice(L);
    lattice.fill(-1.0);

    DEBUG(lattice.hopping_m);

    configuration_t config(lattice, 1.0, U, mu, mu+e_f);
    configuration_t config2(config);
    configuration_t config3(configuration_t(lattice, 1.0, U, mu, mu+e_f));
    config3 = config2;
    config3 = configuration_t(lattice, 1.0, U, mu, mu+e_f);

    Eigen::ArrayXi my_config(lattice.get_msize()); my_config.setZero();
    for (size_t x=0; x<L; x+=1)
        for (size_t y=0; y<L; y+=2) {
            my_config(lattice.pos_to_index({x,y+x%2}))=1;
        }
    config.f_config_ = my_config;
    DEBUG(config.f_config_.transpose());
    config.calc_hamiltonian();
    config.calc_ed();
    DEBUG(config.ed_data().cached_spectrum.transpose());
    //DEBUG(config.cached_evecs);

   }
  catch(triqs::runtime_error const & e) { std::cout  << "exception "<< e.what() << std::endl;}
  return 0;

    return EXIT_SUCCESS;
}
