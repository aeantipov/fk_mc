#ifndef __FK_MC_CONFIGURATION_HPP_
#define __FK_MC_CONFIGURATION_HPP_

#include "common.hpp"
#include <triqs/mc_tools/random_generator.hpp>
#include <numeric>
//#include <triqs/arrays/linalg/eigenelements.hpp>

#include <Eigen/Eigenvalues>

namespace fk {

template <class lattice_t>
struct configuration {

    typedef typename lattice_t::sparse_m sparse_m;
    typedef Eigen::MatrixXd dense_m;
    typedef Eigen::ArrayXi int_array_t;
    typedef Eigen::ArrayXd real_array_t;

    enum class calc_eval { full, arpack };

    const lattice_t& lattice;
    int_array_t f_config;
    double beta, U, mu_c, mu_f;

    mutable sparse_m hamilt;
    size_t n_calc_evals; // number of eigenvalues to calculate
    double eval_weight_tolerance = std::numeric_limits<double>::epsilon();
    mutable real_array_t cached_spectrum;
    mutable real_array_t cached_weights;

    configuration(
        const lattice_t &lattice, double beta, double U, double mu_c, double mu_f):
            lattice(lattice),
            f_config(lattice.m_size),
            beta(beta),
            U(U),mu_c(mu_c),mu_f(mu_f)
            { f_config.setZero(); n_calc_evals = lattice.m_size; }

    void swap(configuration &rhs) {
        f_config.swap(rhs.f_config); 
        cached_spectrum.swap(rhs.cached_spectrum);
        cached_weights.swap(rhs.cached_weights);
        hamilt.swap(rhs.hamilt);
        n_calc_evals = rhs.n_calc_evals;
        beta = rhs.beta; U = rhs.U, mu_c = rhs.mu_c; mu_f = rhs.mu_f; 
        eval_weight_tolerance = rhs.eval_weight_tolerance;
        };

    configuration(const configuration& rhs) = default;
    configuration& operator=(const configuration& rhs) {
        f_config = rhs.f_config; 
        cached_spectrum = rhs.cached_spectrum;
        cached_weights = rhs.cached_weights;
        hamilt = rhs.hamilt;
        n_calc_evals = rhs.n_calc_evals;
        beta = rhs.beta; U = rhs.U, mu_c = rhs.mu_c; mu_f = rhs.mu_f; 
        eval_weight_tolerance = rhs.eval_weight_tolerance;
        return *this;
        };
    configuration(configuration&& rhs):lattice(rhs.lattice) { rhs.swap(*this); };
    configuration& operator=(configuration&& rhs) { rhs.swap(*this); return *this; };

    size_t get_nf() const;
    void randomize_f(triqs::mc_tools::random_generator &rnd, size_t nf = 0);
    size_t get_m_size() const { return lattice.m_size; };

    sparse_m calc_hamiltonian();

    real_array_t  calc_spectrum(calc_eval flag = calc_eval::full);
    real_array_t  calc_spectrum_full();
    //real_array_t  calc_spectrum_arpack();
};

template <class lattice_t>
inline size_t configuration<lattice_t>::get_nf() const
{
    return std::accumulate(f_config.data(), f_config.data()+lattice.m_size, 0);
}

template <class lattice_t>
void configuration<lattice_t>::randomize_f(triqs::mc_tools::random_generator &rnd, size_t nf){
    if (!nf) nf = rnd(lattice.m_size);
    f_config.setZero();
    for (size_t i=0; i<nf; ++i) {  
    size_t ind = rnd(lattice.m_size);
    while (f_config(ind)==1) ind = rnd(lattice.m_size);
    f_config(ind) = 1; 
    };
}

template <class lattice_t>
inline typename configuration<lattice_t>::sparse_m configuration<lattice_t>::calc_hamiltonian()
{
    hamilt.reserve(lattice.nonzero_elems + lattice.m_size);
    hamilt = lattice.hopping_m;
    for (size_t i=0; i<lattice.m_size; ++i) hamilt.coeffRef(i,i)+= -mu_c + U*f_config(i); // unoptimized
    return hamilt;
}


template <class lattice_t>
inline typename configuration<lattice_t>::real_array_t configuration<lattice_t>::calc_spectrum(calc_eval flag)
{
    if (n_calc_evals == lattice.m_size) flag = calc_eval::full; 
    switch (flag) {
        //case calc_eval::arpack: calc_spectrum_arpack(); break;
        case calc_eval::full: calc_spectrum_full(); break;
    }

    cached_weights.resize(cached_spectrum.size());
    n_calc_evals = 0;
    double e0 = cached_spectrum[0]; // unsafe - .minCoeff() is safer but slower;
    double weight0 = exp(-beta*e0);
    //MY_DEBUG("w0 : " << weight0);
    for (size_t i=0; i<cached_spectrum.size(); ++i) { 
        cached_weights(i) = exp(-beta*(cached_spectrum(i)-e0)); 
        if (cached_weights(i)*weight0 >= eval_weight_tolerance) n_calc_evals++; 
        //MY_DEBUG(cached_weights(i) << " " << exp(-beta*(cached_spectrum(i))) <<  " " << eval_weight_tolerance << " " << n_calc_evals);
        };
    n_calc_evals++; // calculate extra eigenvalue below the tolerance for future checks
    if (n_calc_evals > lattice.m_size) 
        n_calc_evals = lattice.m_size;
        else if (n_calc_evals == cached_spectrum.size()+1) 
            n_calc_evals++;
    //MY_DEBUG(n_calc_evals << " " << cached_spectrum.size() << "|" << eval_weight_tolerance);
    return cached_spectrum;
}

template <class lattice_t>
inline typename configuration<lattice_t>::real_array_t configuration<lattice_t>::calc_spectrum_full()
{
    dense_m h(hamilt);
    Eigen::SelfAdjointEigenSolver<dense_m> s(h,Eigen::EigenvaluesOnly);
    cached_spectrum = s.eigenvalues();
    return cached_spectrum;
}

/* Obsolete, clean 
template <class lattice_t>
inline typename configuration<lattice_t>::real_array_t configuration<lattice_t>::calc_spectrum_arpack()
{

    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<sparse_m> s(hamilt,n_calc_evals,"SA",Eigen::EigenvaluesOnly);
    cached_spectrum = s.eigenvalues().reverse();
    return cached_spectrum;
}
*/



} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_CONFIGURATION_HPP_

