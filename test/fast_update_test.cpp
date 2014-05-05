#include <boost/mpi/environment.hpp>
#include <chrono>
#include <Eigen/Sparse>

#include "eigen/ArpackSupport"
//#include <triqs/gfs.hpp>

#include "fk_mc.hpp"
#include "moves.hpp"
#include "chebyshev.hpp"

using namespace fk;

typedef configuration_t::dense_m dense_m;
typedef configuration_t::sparse_m sparse_m;

//typedef triqs::gfs::linear_mesh<triqs::gfs::R_domain> grid_type;
//grid_type cr_grid(triqs::gfs::R_domain(), -1+0.001, 1-0.001, 1000, triqs::gfs::full_bins);

int main(int argc, char* argv[])
{
    boost::mpi::environment env(argc, argv);
  try {

    size_t L = 16;

    double U = 1.5;
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
    config.calc_ed();

    configuration_t config1(config);
    config1 = config;

    move_addremove move1(beta, config, r);
    auto w = move1.attempt();
    move1.accept();

    /*std::cout << "Old config  : " << config1.f_config_.transpose() << std::endl;
    std::cout << "New config  : " << config.f_config_.transpose() << std::endl;
    std::cout << "Old spectrum: " << config1.cached_spectrum.transpose() << std::endl;
    std::cout << "New spectrum: " << config.cached_spectrum.transpose() + 0.0*(config1.cached_spectrum[0] - config.cached_spectrum[0]) << std::endl;
    */

    sparse_m h = config1.hamilt_;
    auto e_min = Eigen::ArpackGeneralizedSelfAdjointEigenSolver<sparse_m>(h,1,"SA",Eigen::EigenvaluesOnly).eigenvalues()[0];
    auto e_max = Eigen::ArpackGeneralizedSelfAdjointEigenSolver<sparse_m>(h,1,"LA",Eigen::EigenvaluesOnly).eigenvalues()[0];
    std::cout << std::setprecision(5) <<  std::setw(3) << "h : lowest eigenvalue :" << e_min << " top eigenvalue: " << e_max << std::endl;

    double a = (e_max - e_min)/2.;
    double b = (e_max + e_min)/2.; 

    MY_DEBUG((a*(-1.) + b) << " == " << e_min);
    MY_DEBUG((a*(1.) + b) << " == " << e_max);


    sparse_m x = h; 
    for (size_t i=0; i<lattice.get_msize(); ++i) x.coeffRef(i,i)+= -b; // unoptimized
    x/=a;
    e_min = Eigen::ArpackGeneralizedSelfAdjointEigenSolver<sparse_m>(x,1,"SA",Eigen::EigenvaluesOnly).eigenvalues()[0];
    e_max = Eigen::ArpackGeneralizedSelfAdjointEigenSolver<sparse_m>(x,1,"LA",Eigen::EigenvaluesOnly).eigenvalues()[0];
    std::cout << std::setprecision(5) <<  std::setw(3) << "x : lowest eigenvalue :" << e_min << " top eigenvalue: " << e_max << std::endl;

    dense_m basis = dense_m::Identity(lattice.get_msize(),lattice.get_msize());

    size_t cheb_size = int(std::log(lattice.get_msize())*1.5);

    std::vector<dense_m> cm(cheb_size);
    std::vector<double> moments(cheb_size, 0);
    cm[0] = basis;
    cm[1] = x*basis;

    //MY_DEBUG(cr_grid[0]);
    //MY_DEBUG(cr_grid[500]);
    //MY_DEBUG(cr_grid.nearest_index(0.0));

    auto s_f = [a,b,beta,&lattice](double w){return lattice.get_msize()*log(1. + exp(-beta*(a*w+b)));}; 
    auto x_f = [](double x){return x;};
    auto x2_f = [](double x){return std::pow(sin(-x*M_PI),2);};

    Eigen::VectorXd chebyshev_grid = Eigen::VectorXd::LinSpaced(Eigen::Sequential, 2000, -1.+0.001, 1-0.001) ;
    size_t ngrid_points = 2000;
    chebyshev::chebyshev_eval ch(cheb_size, ngrid_points);
    //for (size_t i=0; i<cheb_grid2.size(); ++i) cheb_grid2[i] = -cos(M_PI * (2.*i + 1.) / (2.*ngrid_points));
    //MY_DEBUG(std::setprecision(12) << cheb_grid2.transpose());


    MY_DEBUG(chebyshev_grid.head(1) << " " << chebyshev_grid.tail(1));

    double s=0;
    for (int i=0; i<cheb_size; i++) {
            //sparse_m t = chebyshev_t(x,i);
            //dense_m vm = t*basis;
            if (i>1) cm[i] = x*2*cm[i-1] - cm[i-2];
            moments[i] = cm[i].diagonal().sum()/lattice.get_msize();
            std::cout << "moment [" << i << "] = " << moments[i] << std::endl;
            auto s_m2 = ch.moment(s_f, i);
            s+=(i==0?1.:2.)*s_m2*moments[i];
            std::cout << s_m2 << " | " << s << std::endl;
        }

    MY_DEBUG(s);
    double s2 = 0.0;
    for (int i=0; i<lattice.get_msize(); i++) { double e = config1.ed_data().cached_spectrum[i]; s2 += std::log(1 + exp(-beta * e)); };
    MY_DEBUG(s2);
       }
  catch(triqs::runtime_error const & e) { std::cout  << "exception "<< e.what() << std::endl; return EXIT_FAILURE;}

    return EXIT_SUCCESS;
}
