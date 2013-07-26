#include "fk_mc.hpp"
#include <boost/mpi/environment.hpp>
#include <chrono>

using namespace fk;

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
  try {
    size_t L = 4;

    double U = 6.0;
    double t = -1.0;
    double T = 0.15;
    double beta = 1.0/T;

    square_lattice_traits<2> lattice(L);
    lattice.fill(-1.0);

    DEBUG(lattice.get_hopping_matrix());

    fk_mc<square_lattice_traits<2>> mc(lattice);

    triqs::utility::parameters p;
    p["U"] = U;
    p["beta"] = beta;
    p["Nf_start"] = L*L/2;
    p["Random_Generator_Name"] = ""; 
    p["Random_Seed"] = 34788+std::chrono::system_clock::now().time_since_epoch().count();
    p["Verbosity"] = 3;
    p["Length_Cycle"] = 15; 
    p["N_Warmup_Cycles"] = 1;
    p["N_Cycles"] = 1;
    p["max_time"]=5;

    mc.solve(p);
  }
  catch(triqs::runtime_error const & e) { std::cout  << "exception "<< e.what() << std::endl;}
  return 0;

    return EXIT_SUCCESS;
}
