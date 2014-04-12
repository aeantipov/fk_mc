#include <boost/mpi/environment.hpp>
#include <chrono>
#include <Eigen/Sparse>

#include "eigen/ArpackSupport"
#include <triqs/gfs.hpp>

#include "fk_mc.hpp"
#include "moves.hpp"


using namespace fk;

typedef configuration_t::dense_m dense_m;
typedef configuration_t::sparse_m sparse_m;

double pi = std::acos(-1);
//typedef triqs::gfs::linear_mesh<triqs::gfs::R_domain> grid_type;
//grid_type cr_grid(triqs::gfs::R_domain(), -1+0.001, 1-0.001, 1000, triqs::gfs::full_bins);

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
  try {

    size_t L = 24;

    double U = 1.0;
    double mu = U/2;
    double e_f = 0.0;
    double t = -1.0;
    double T = 0.16;
    double beta = 1.0/T;

    triqs::mc_tools::random_generator r("mt19937", 32167);

    typedef hypercubic_lattice<2> lattice_t;
    lattice_t lattice(L);
    lattice.fill(-1.0);

    configuration_t config(lattice, 1.0, U, mu, mu+e_f);
    config.randomize_f(r,L*L/2);
    config.calc_hamiltonian();
    config.calc_full_spectrum(true);

    configuration_t config1(config);
    config1 = config;

    move_addremove move1(beta, config, r);
    auto w = move1.attempt();
    move1.accept();

    /*std::cout << "Old config  : " << config1.f_config.transpose() << std::endl;
    std::cout << "New config  : " << config.f_config.transpose() << std::endl;
    std::cout << "Old spectrum: " << config1.cached_spectrum.transpose() << std::endl;
    std::cout << "New spectrum: " << config.cached_spectrum.transpose() + 0.0*(config1.cached_spectrum[0] - config.cached_spectrum[0]) << std::endl;
    */

    sparse_m h = config1.hamilt;
    auto e_min = Eigen::ArpackGeneralizedSelfAdjointEigenSolver<sparse_m>(h,1,"SA",Eigen::EigenvaluesOnly).eigenvalues()[0];
    auto e_max = Eigen::ArpackGeneralizedSelfAdjointEigenSolver<sparse_m>(h,1,"LA",Eigen::EigenvaluesOnly).eigenvalues()[0];
    std::cout << std::setprecision(5) <<  std::setw(3) << "h : lowest eigenvalue :" << e_min << " top eigenvalue: " << e_max << std::endl;

    double a = (e_max - e_min)/2.;
    double b = (e_max + e_min)/2.; 

    sparse_m x = h; 
    for (size_t i=0; i<lattice.get_msize(); ++i) x.coeffRef(i,i)+= -b; // unoptimized
    x/=a;
    e_min = Eigen::ArpackGeneralizedSelfAdjointEigenSolver<sparse_m>(x,1,"SA",Eigen::EigenvaluesOnly).eigenvalues()[0];
    e_max = Eigen::ArpackGeneralizedSelfAdjointEigenSolver<sparse_m>(x,1,"LA",Eigen::EigenvaluesOnly).eigenvalues()[0];
    std::cout << std::setprecision(5) <<  std::setw(3) << "x : lowest eigenvalue :" << e_min << " top eigenvalue: " << e_max << std::endl;

    dense_m basis = dense_m::Identity(lattice.get_msize(),lattice.get_msize());

    size_t cheb_size = 200;

    std::vector<dense_m> cm(cheb_size);
    std::vector<double> moments(cheb_size, 0);
    cm[0] = basis;
    cm[1] = x*basis;

    //MY_DEBUG(cr_grid[0]);
    //MY_DEBUG(cr_grid[500]);
    //MY_DEBUG(cr_grid.nearest_index(0.0));

    auto s_f = [a,b,beta](double w){return log(1. + exp(-beta*(a*w+b)));}; 
    auto x_f = [](double x){return x;};
    auto x2_f = [](double x){return x*x;};

    Eigen::VectorXd chebyshev_grid = Eigen::VectorXd::LinSpaced(Eigen::Sequential, 1000, -1.+0.001, 1-0.001) ;

    MY_DEBUG(simpson(x_f,chebyshev_grid));
    MY_DEBUG(simpson(x2_f, chebyshev_grid));

    double s=1;
    for (int i=0; i<cheb_size; i++) {
            //sparse_m t = chebyshev_t(x,i);
            //dense_m vm = t*basis;
            if (i>1) cm[i] = x*2*cm[i-1] - cm[i-2];
            moments[i] = cm[i].diagonal().sum()/lattice.get_msize();
            std::cout << "moment [" << i << "] = " << moments[i] << std::endl;
            auto s_m = chebyshev_moment(s_f, i, chebyshev_grid);
            s+=2.*s_m;
            std::cout << s_m << std::endl;
        }

    MY_DEBUG(s);

    if (0==1) { 
        dense_m a = config.hamilt;
        dense_m diag(a); a.setZero(); a.diagonal().setOnes();
        MY_DEBUG((a + diag*2.) * (a - diag*2.).inverse())
        }

    if (0==1) { 
    dense_m h0 = config.hamilt;
    dense_m v = config1.hamilt - config.hamilt;

    auto commutator = [](dense_m x, dense_m y)->dense_m { return x*y - y*x; };
    MY_DEBUG("V = " << v);
    MY_DEBUG("[H,V] = " << commutator (h0, v));
    MY_DEBUG("[H,[H,V]] = " << commutator(h0,commutator (h0, v)));
    MY_DEBUG("[V,[H,V]] = " << commutator(v,commutator (h0, v)));
    };
   }
  catch(triqs::runtime_error const & e) { std::cout  << "exception "<< e.what() << std::endl; return EXIT_FAILURE;}

    return EXIT_SUCCESS;
}
