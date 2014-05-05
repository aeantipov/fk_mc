#ifndef __FK_MC_CONFIGURATION_HPP_
#define __FK_MC_CONFIGURATION_HPP_

#include <triqs/mc_tools/random_generator.hpp>
#include <numeric>

#include <Eigen/Eigenvalues>

#include "common.hpp"
#include "lattice.hpp"
#include "chebyshev.hpp"

namespace fk {

struct config_params {
    double beta, U, mu_c, mu_f;
    bool operator== ( const config_params& rhs) const;
};

struct ed_cache { 
    enum status_eval {empty, spectrum, full};
    typedef typename lattice_base::sparse_m sparse_m;
    typedef Eigen::ArrayXd real_array_t;
    typedef Eigen::MatrixXd dense_m;
    
    status_eval status;
    real_array_t cached_spectrum;
    dense_m cached_evecs;
    real_array_t cached_weights;
};

struct chebyshev_cache { 
    typedef typename lattice_base::sparse_m sparse_m;
    enum status_eval {empty, full};

    status_eval status;
    double e_max;
    double e_min;
    double a; // (e_max - e_min)/2.
    double b; // (e_max + e_min)/2.
    sparse_m x;
    std::vector<double> moments;
};

struct configuration_t {
    typedef typename ed_cache::sparse_m sparse_m;
    typedef typename ed_cache::dense_m dense_m;
    typedef typename ed_cache::real_array_t real_array_t;
    typedef Eigen::ArrayXi int_array_t;

    configuration_t(const lattice_base &lattice, double beta, double U, double mu_c, double mu_f);

    configuration_t(const configuration_t& rhs) = default ;
    configuration_t& operator=(const configuration_t& rhs);
    configuration_t(configuration_t&& rhs) = default;
    configuration_t& operator=(configuration_t&& rhs) = default;
    void swap(configuration_t &rhs);

    size_t get_nf() const;
    void randomize_f(triqs::mc_tools::random_generator &rnd, size_t nf = 0);
    const sparse_m& calc_hamiltonian();
    void reset_cache(){ed_data_.status =  ed_cache::empty; cheb_data_.status = chebyshev_cache::empty;}

    real_array_t calc_ed(bool calc_evecs = false);
    void calc_chebyshev(const chebyshev_eval& cheb);

    const config_params& params(){return params_;}
///
    const lattice_base& lattice_;
    const config_params params_;
    int_array_t f_config_;
    sparse_m hamilt_;
    ed_cache ed_data_;
    chebyshev_cache cheb_data_;
};

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_CONFIGURATION_HPP_

