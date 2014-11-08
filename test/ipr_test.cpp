#include <boost/mpi/environment.hpp>
#include <chrono>
#include "lattice/hypercubic.hpp"
#include "fk_mc.hpp"
#include "measures/ipr.hpp"

using namespace fk;

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator comm;
  try {

    size_t L = 4;

    double U = 1.0;
    double mu = U/2;
    double e_f = 0.0;
    double t = -1.0;
    double T = 0.15;
    double beta = 1.0/T;

    typedef hypercubic_lattice<2> lattice_t;
    lattice_t lattice(L);
    lattice.fill(-1.0);

    FKDEBUG(lattice.hopping_m());

    configuration_t config(lattice, 1.0, U, mu, mu+e_f);
    configuration_t config2(config);
    configuration_t config3(configuration_t(lattice, 1.0, U, mu, mu+e_f));
    config3 = config2;
    config3 = configuration_t(lattice, 1.0, U, mu, mu+e_f);

    Eigen::ArrayXi my_config(lattice.get_msize()); my_config.setZero();
    for (int x=0; x<L; x+=1)
        for (int y=0; y<L; y+=2) {
            my_config(lattice.pos_to_index({{x,y+x%2}}))=1;
        }
    config.f_config_ = my_config;
    FKDEBUG(config.f_config_.transpose());
    config.calc_hamiltonian();
    config.calc_ed(true);
    FKDEBUG(config.ed_data().cached_spectrum.transpose());
    //FKDEBUG(config.cached_evecs);


    std::vector<std::vector<double>> ipr_vals;
    auto ipr_m = measure_ipr<lattice_t>(config,lattice,ipr_vals);

    ipr_m.accumulate(1.0);
    ipr_m.accumulate(1.0);
    ipr_m.collect_results(comm);
    for (auto& v:ipr_vals) std::cout << v[1] << " " << std::flush; std::cout << std::endl;

   }
  catch(triqs::runtime_error const & e) { std::cout  << "exception "<< e.what() << std::endl; return EXIT_FAILURE;}

  return EXIT_SUCCESS;
}
