#include <boost/mpi/environment.hpp>
#include <chrono>
#include <Eigen/Sparse>

#include <chrono>
#include<random>
 
using namespace std::chrono;
//using std::chrono::duration_cast;
//using std::chrono::microseconds;
//using std::chrono::steady_clock;

#include <gtest/gtest.h>

#include "eigen/ArpackSupport"

#include "lattice/hypercubic.hpp"
#include "fk_mc.hpp"
#include "moves.hpp"
#include "chebyshev.hpp"
#include "moves_chebyshev.hpp"

#include <tclap/CmdLine.h>

using namespace fk;

typedef configuration_t::dense_m dense_m;
typedef configuration_t::sparse_m sparse_m;

double weight_tol = 6e-2;

int L;
double U,T,cheb_prefactor;
size_t rnd_seed;

TEST(FastUpdateTest, weight) { 
    std::cout << "T = " << T <<  "; U = " << U << "; L = " << L << std::endl;
    std::cout << "cheb prefactor : " << cheb_prefactor << std::endl;

    double mu = U/2.;
    double e_f = 0.0;
    double t = 1.0;
    double beta = 1.0/T;

    typedef hypercubic_lattice<2> lattice_t;
    lattice_t lattice(L);
    lattice.fill(-1.0);

    std::cout << "Random seed : " << rnd_seed << std::endl;
    std::cout << "Number of states = " << lattice.msize() << std::endl;
    std::cout << "log(nstates) = " << int(std::log(lattice.msize())) << std::endl;

    size_t cheb_size = std::min( int(std::log(lattice.msize())*cheb_prefactor), lattice.msize()/4);
    cheb_size+=cheb_size%2;
    std::cout << "Number of Chebyshev polynomials = " << cheb_size << std::endl;
    size_t ngrid_points = cheb_size*2;// std::max(size_t(200), cheb_size);
    chebyshev::chebyshev_eval ch(cheb_size, ngrid_points);

    configuration_t config(lattice, beta, U, mu, mu+e_f);
    //triqs::mc_tools::random_generator r1("mt19937", rnd_seed);
    random_generator r1 (rnd_seed);
    config.randomize_f(r1,L*L/2);
    config.calc_hamiltonian();
    config.calc_ed();
    config.calc_chebyshev(ch);

    auto correct_w = [](double w1, double w2){ 
        return 
        (std::abs(w1) + std::abs(w1) < 2.*weight_tol) ||
        (std::abs((w1 - w2))<weight_tol) ||
        (std::abs(w1) > 1 - weight_tol && std::abs(w2) > 1-weight_tol);
        };

   // test move
    //triqs::mc_tools::random_generator r("mt19937", rnd_seed);
    random_generator r(rnd_seed);
    random_generator r2(r);
    //triqs::mc_tools::random_generator r2(r);
    bool result;
    steady_clock::time_point start, end;

    configuration_t config1(config), config2(config);
    move_addremove move1(beta, config1, r);
    start = steady_clock::now();
    auto w = move1.attempt();
    end = steady_clock::now();
    std::cout << "move add_remove (ed)   weight : " << w << std::endl;
    std::cout << "time duration : " << duration_cast<milliseconds>(end-start).count() << " ms" << std::endl;

    chebyshev::move_addremove move_c(beta, config2, ch, r2); 
    start = steady_clock::now();
    auto w_c = move_c.attempt();
    end = steady_clock::now();
    std::cout << "move add_remove (cheb) weight : " << w_c << std::endl;
    std::cout << "time duration : " << duration_cast<milliseconds>(end-start).count() << " ms" << std::endl;

    ASSERT_EQ(((config2.f_config_ - config1.f_config_).sum()),0);
    result = correct_w(w,w_c);
    if (result) std::cout << "-->weight diff: " << std::abs((w_c - w)) << " ; tol = " << weight_tol << std::endl;
    EXPECT_EQ(result, true);


    start = steady_clock::now();
    w = move_flip(beta, config1, r).attempt();
    end = steady_clock::now();
    std::cout << "move flip (ed)   weight : " << w << std::endl;
    std::cout << "time duration : " << duration_cast<milliseconds>(end-start).count() << " ms" << std::endl;

    start = steady_clock::now();
    w_c = chebyshev::move_flip(beta, config2, ch, r2).attempt();
    end = steady_clock::now();
    std::cout << "move flip (cheb) weight : " << w_c << std::endl;
    std::cout << "time duration : " << duration_cast<milliseconds>(end-start).count() << " ms" << std::endl;

    ASSERT_EQ(((config2.f_config_ - config1.f_config_).sum()),0);
    result = correct_w(w,w_c);
    if (result) std::cout << "-->weight diff: " << std::abs((w_c - w)) << " ; tol = " << weight_tol << std::endl;
    EXPECT_EQ(result, true);


    start = steady_clock::now();
    w = move_randomize(beta, config1, r).attempt();
    end = steady_clock::now();
    std::cout << "move randomize (ed)   weight : " << w << std::endl;
    std::cout << "time duration : " << duration_cast<milliseconds>(end-start).count() << " ms" << std::endl;

    start = steady_clock::now();
    w_c = chebyshev::move_randomize(beta, config2, ch, r2).attempt();
    end = steady_clock::now();
    std::cout << "move randomize (cheb) weight : " << w_c << std::endl;
    std::cout << "time duration : " << duration_cast<milliseconds>(end-start).count() << " ms" << std::endl;

    ASSERT_EQ(((config2.f_config_ - config1.f_config_).sum()),0);
    result = correct_w(w,w_c);
    if (result) std::cout << "-->weight diff: " << std::abs((w_c - w)) << " ; tol = " << weight_tol << std::endl;
    EXPECT_EQ(result, true);
    std::cout << "Random seed : " << rnd_seed << std::endl;
}


int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    boost::mpi::communicator world;

    TCLAP::CmdLine cmd("Falicov-Kimball Monte Carlo weight check", ' ', "");
    TCLAP::ValueArg<double> U_arg("U","U","value of U",false,1.0,"double",cmd);
    TCLAP::ValueArg<double> T_arg("T","T","Temperature",false,0.1,"double",cmd);
    TCLAP::ValueArg<size_t> L_arg("L","L","system size",false,16,"int",cmd);
    TCLAP::ValueArg<double>   cheb_pref_arg("c","cheb_prefactor","Prefactor of log(N) chebyshev polynomials", false, 2.35, "double", cmd);
    TCLAP::SwitchArg          fixed_seed_switch("f","f","Make a random or fixed seed?", cmd, false);
    TCLAP::ValueArg<int>   seed_arg("","seed","Random seed (otherwise random)", false, 32167, "int", cmd);
    cmd.parse( argc, argv );

    U = U_arg.getValue();
    T = T_arg.getValue();
    L = L_arg.getValue();
    cheb_prefactor = cheb_pref_arg.getValue();
    rnd_seed = (!fixed_seed_switch.getValue()?std::random_device()():(seed_arg.getValue()+world.rank())); 

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
