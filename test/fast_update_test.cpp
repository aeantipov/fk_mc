#include <boost/mpi/environment.hpp>
#include <chrono>
#include <Eigen/Sparse>

#include <chrono>
 
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

double weight_tol = 4e-2;

TEST(Chebyshev,diff)
{
    size_t L = 16;

    double U = 8;
    double mu = U/2;
    double e_f = 0.0;
    double t = -1.0;
    double T = 0.16;
    double beta = 1.0/T;


    typedef hypercubic_lattice<2> lattice_t;
    lattice_t lattice(L);
    lattice.fill(-1.0);

    configuration_t config(lattice, beta, U, mu, mu+e_f);
    triqs::mc_tools::random_generator r1("mt19937", 32167);
    config.randomize_f(r1,L*L/2);
    config.calc_hamiltonian();
    config.calc_ed();

    std::cout << "Number of states = " << lattice.get_msize() << std::endl;
    std::cout << "log(nstates) = " << int(std::log(lattice.get_msize())) << std::endl;
    size_t cheb_size = int(std::log(lattice.get_msize()))*2;
    cheb_size+=cheb_size%2;
    size_t ngrid_points = cheb_size; //std::max(size_t(200), cheb_size);
    chebyshev::chebyshev_eval ch(cheb_size, ngrid_points);
    cheb_size = ch.cheb_size();
    std::cout << "Number of Chebyshev polynomials = " << cheb_size << std::endl;

    ASSERT_NEAR(ch.moment([](double x){return 1.;},0),1., 1e-8);
    ASSERT_NEAR(ch.moment([](double x){return x*x;},0),0.5, 1e-8);
    ASSERT_NEAR(ch.moment([](double x){return sin(x);},1),0.440051, 1e-6);
    ASSERT_NEAR(ch.moment([](double x){return sin(x);},5),0.000249758, 1e-9);

    config.calc_chebyshev(ch);
    auto logZ = config.cheb_data_.logZ;
    
    config.calc_ed();
    double s2 = 0.0;
    for (int i=0; i<lattice.get_msize(); i++) { double e = config.ed_data().cached_spectrum[i]; s2 += std::log(1 + exp(-beta * e)); };

    std::cout << "logZ chebyshev : " << logZ << std::endl;
    std::cout << "logZ ed        : " << s2 << std::endl;
    std::cout << "logZ ed        : " << config.ed_data_.logZ << std::endl;
    EXPECT_EQ( (std::abs((logZ - s2) / logZ) > 1e-2), false);
    EXPECT_EQ( (std::abs((config.ed_data_.logZ - s2) / logZ) > 1e-15), false);
}


int main(int argc, char* argv[])
{

    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
