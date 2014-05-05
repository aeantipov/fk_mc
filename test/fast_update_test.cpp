#include <boost/mpi/environment.hpp>
#include <chrono>
#include <Eigen/Sparse>

#include "eigen/ArpackSupport"
//#include <triqs/gfs.hpp>

#include "fk_mc.hpp"
#include "moves.hpp"
#include "chebyshev.hpp"
#include "moves_chebyshev.hpp"

using namespace fk;

typedef configuration_t::dense_m dense_m;
typedef configuration_t::sparse_m sparse_m;

//typedef triqs::gfs::linear_mesh<triqs::gfs::R_domain> grid_type;
//grid_type cr_grid(triqs::gfs::R_domain(), -1+0.001, 1-0.001, 1000, triqs::gfs::full_bins);

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
    size_t L = 20;

    double U = 1.5;
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

    size_t cheb_size = int(std::log(lattice.get_msize())*1.5);
    size_t ngrid_points = std::max(size_t(200), cheb_size);
    chebyshev::chebyshev_eval ch(cheb_size, ngrid_points);

    config.calc_chebyshev(ch);
    auto logZ = config.cheb_data_.logZ;
    
    config.calc_ed();
    double s2 = 0.0;
    for (int i=0; i<lattice.get_msize(); i++) { double e = config.ed_data().cached_spectrum[i]; s2 += std::log(1 + exp(-beta * e)); };

    std::cout << "logZ chebyshev : " << logZ << std::endl;
    std::cout << "logZ ed        : " << s2 << std::endl;
    std::cout << "logZ ed        : " << config.ed_data_.logZ << std::endl;
    if (std::abs((logZ - s2) / logZ) > 1e-2) return EXIT_FAILURE;
    if (std::abs((config.ed_data_.logZ - s2) / logZ) > 1e-15) return EXIT_FAILURE;

    // test move
    triqs::mc_tools::random_generator r("mt19937", 32167);
    triqs::mc_tools::random_generator r2(r);

    configuration_t config1(config), config2(config);
    move_addremove move1(beta, config1, r);
    auto w = move1.attempt();

    chebyshev::move_addremove move_c(beta, config2, ch, r2); 
    auto w_c = move_c.attempt();

    if ((config2.f_config_ - config1.f_config_).sum()!=0) return EXIT_FAILURE;

    std::cout << "move add_remove (ed)   weight : " << w << std::endl;
    std::cout << "move add_remove (cheb) weight : " << w_c << std::endl;

    w = move_flip(beta, config1, r).attempt();
    w_c = chebyshev::move_flip(beta, config2, ch, r2).attempt();
    if ((config2.f_config_ - config1.f_config_).sum()!=0) return EXIT_FAILURE;
    std::cout << "move flip (ed)   weight : " << w << std::endl;
    std::cout << "move flip (cheb) weight : " << w_c << std::endl;

    w = move_randomize(beta, config1, r).attempt();
    w_c = chebyshev::move_randomize(beta, config2, ch, r2).attempt();
    if ((config2.f_config_ - config1.f_config_).sum()!=0) return EXIT_FAILURE;
    std::cout << "move randomize (ed)   weight : " << w << std::endl;
    std::cout << "move randomize (cheb) weight : " << w_c << std::endl;

    return EXIT_SUCCESS;
}
