#ifndef __FK_MC_CONFIGURATION_HPP_
#define __FK_MC_CONFIGURATION_HPP_

#include <triqs/mc_tools/random_generator.hpp>
#include <numeric>

#include <Eigen/Eigenvalues>

#include "common.hpp"
#include "lattice.hpp"

namespace fk {

struct configuration_t {
    typedef typename lattice_base::sparse_m sparse_m;
    typedef Eigen::MatrixXd dense_m;
    typedef Eigen::ArrayXi int_array_t;
    typedef Eigen::ArrayXd real_array_t;

    const lattice_base& lattice;
    int_array_t f_config;
    double beta, U, mu_c, mu_f;

    mutable sparse_m hamilt;
    double eval_weight_tolerance = std::numeric_limits<double>::epsilon();
    mutable real_array_t cached_spectrum;
    mutable dense_m cached_evecs;
    mutable real_array_t cached_weights;

    configuration_t(
        const lattice_base &lattice, double beta, double U, double mu_c, double mu_f):
            lattice(lattice),
            f_config(lattice.get_msize()),
            hamilt(lattice.hopping_m.rows(), lattice.hopping_m.cols()),
            beta(beta),
            U(U),mu_c(mu_c),mu_f(mu_f)
            { f_config.setZero(); }

    void swap(configuration_t &rhs) {
        f_config.swap(rhs.f_config); 
        cached_spectrum.swap(rhs.cached_spectrum);
        cached_weights.swap(rhs.cached_weights);
        hamilt.swap(rhs.hamilt);
        beta = rhs.beta; U = rhs.U, mu_c = rhs.mu_c; mu_f = rhs.mu_f; 
        eval_weight_tolerance = rhs.eval_weight_tolerance;
        };
    configuration_t(const configuration_t& rhs) = default ;
    configuration_t& operator=(const configuration_t& rhs) {
        f_config = rhs.f_config; 
        cached_spectrum = rhs.cached_spectrum;
        cached_weights = rhs.cached_weights;
        hamilt = rhs.hamilt;
        beta = rhs.beta; U = rhs.U, mu_c = rhs.mu_c; mu_f = rhs.mu_f; 
        eval_weight_tolerance = rhs.eval_weight_tolerance;
        return *this;
        };
    configuration_t(configuration_t&& rhs):lattice(rhs.lattice) { rhs.swap(*this); };
    configuration_t& operator=(configuration_t&& rhs) { rhs.swap(*this); return *this; };
    size_t get_nf() const;
    void randomize_f(triqs::mc_tools::random_generator &rnd, size_t nf = 0);

    sparse_m calc_hamiltonian();

    real_array_t  calc_full_spectrum(bool calc_evecs = false);
};

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_CONFIGURATION_HPP_

