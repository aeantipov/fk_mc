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
    double t = -1.0;
    double T = 0.1;
    double beta = 1.0/T;

    triangular_lattice_traits lattice(L);
    lattice.fill(-1.0,0);

    MY_DEBUG(lattice.hopping_m);

    fk_mc<triangular_lattice_traits > mc(lattice);

    triqs::utility::parameters p;
    p["U"] = U;
    p["mu_c"] = U/2; p["mu_f"] = U/2;
    p["beta"] = beta;
    p["Nf_start"] = L*L/2;
    p["Random_Seed"] = 34788;
    p["Verbosity"] = 3;
    p["Length_Cycle"] = 1; 
    p["N_Warmup_Cycles"] = 1;
    p["N_Cycles"] = 400;
    p["max_time"]=5;
    p["measure_spectrum_history"] = true;

    mc.solve(p);

  }
  catch(triqs::runtime_error const & e) { std::cout  << "exception "<< e.what() << std::endl;}
  return 0;

    return EXIT_SUCCESS;
}
