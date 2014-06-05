#include <boost/mpi/environment.hpp>
#include <chrono>

#include "gtest/gtest.h"

#include "lattice/triangular.hpp"
#include "fk_mc.hpp"
#include "data_saveload.hpp"

using namespace fk;

TEST(fk_mc, save)
{
    boost::mpi::communicator world;
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
    p["random_seed"] = 32167;
    p["verbosity"] = 1;
    p["n_warmup_cycles"] = 1;
    p["n_cycles"] = 10;
    p["length_cycle"] = 2; 
    p["max_time"]=5;
    p["measure_history"] = true;

    fk_mc<triangular_lattice> mc(lattice,p);
    mc.solve();

    world.barrier();
    if (world.rank() == 0) {
        save_data(mc,p,"saveload_test.h5",false);
        }
}

TEST(fk_mc, load)
{
    //fk_mc mc = load_data("saveload_test.h5"); 
}

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
