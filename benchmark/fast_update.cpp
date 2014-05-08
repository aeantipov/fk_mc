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
//#include <triqs/gfs.hpp>

#include "fk_mc.hpp"
#include "moves.hpp"
#include "chebyshev.hpp"
#include "moves_chebyshev.hpp"

using namespace fk;

typedef configuration_t::dense_m dense_m;
typedef configuration_t::sparse_m sparse_m;

double weight_tol = 6e-2;

typedef std::tuple<double,double,int> TUL_tuple;
class FastUpdateTest : public ::testing::TestWithParam<TUL_tuple> {
public:
    size_t seed1 = std::random_device()();
};

std::vector<double> T_vals({{0.15}});
std::vector<double> U_vals({{2.0, 16.0}});
std::vector<int>    L_vals({{32,48}});

std::vector<TUL_tuple> generate() { 
    std::vector<TUL_tuple> out;
    out.reserve(T_vals.size() * U_vals.size() * L_vals.size());
    for (double T:T_vals) 
        for (double U:U_vals) 
            for (int L:L_vals) 
                out.push_back(std::make_tuple(T,U,L));
    return out;
}

TEST_P(FastUpdateTest, weight) { 
    double T,U; int L;
    std::tie(T,U,L) = GetParam();
    std::cout << "T = " << T <<  "; U = " << U << "; L = " << L << std::endl;

    double mu = U/2.;
    double e_f = 0.0;
    double t = -1.0;
    double beta = 1.0/T;

    typedef hypercubic_lattice<2> lattice_t;
    lattice_t lattice(L);
    lattice.fill(-1.0);

    std::cout << "Number of states = " << lattice.get_msize() << std::endl;
    std::cout << "log(nstates) = " << int(std::log(lattice.get_msize())) << std::endl;

    double cheb_prefactor = 2.5 ;

    size_t cheb_size = std::min( int(std::log(lattice.get_msize())*cheb_prefactor), lattice.get_msize()/4);
    cheb_size+=cheb_size%2;
    std::cout << "Number of Chebyshev polynomials = " << cheb_size << std::endl;
    size_t ngrid_points = cheb_size*2;// std::max(size_t(200), cheb_size);
    chebyshev::chebyshev_eval ch(cheb_size, ngrid_points);

    configuration_t config(lattice, beta, U, mu, mu+e_f);
    triqs::mc_tools::random_generator r1("mt19937", seed1);
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
    triqs::mc_tools::random_generator r("mt19937", seed1);
    triqs::mc_tools::random_generator r2(r);
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
}

INSTANTIATE_TEST_CASE_P(ULtest,
                        FastUpdateTest,::testing::ValuesIn(generate()));

int main(int argc, char* argv[])
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
