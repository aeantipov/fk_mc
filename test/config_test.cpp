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

    square_lattice_traits<2> lattice(L);
    lattice.fill(-1.0);

    DEBUG(lattice.get_hopping_matrix());

    configuration<square_lattice_traits<2>> config(lattice, U, mu, mu+e_f);

    real_array_t my_config(lattice.m_size); my_config()=0;
    for (size_t x=0; x<L; x+=1)
        for (size_t y=0; y<L; y+=2) {
            my_config(lattice.pos_to_index({x,y+x%2}))=1;
        }
    config.f_config = my_config;
    DEBUG(config.f_config);
    DEBUG(config.get_spectrum());

   }
  catch(triqs::runtime_error const & e) { std::cout  << "exception "<< e.what() << std::endl;}
  return 0;

    return EXIT_SUCCESS;
}
