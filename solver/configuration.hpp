#ifndef __FK_MC_CONFIGURATION_HPP_
#define __FK_MC_CONFIGURATION_HPP_

#include "common.hpp"
#include <triqs/mc_tools/random_generator.hpp>
#include <numeric>
#include <triqs/arrays/linalg/eigenelements.hpp>

#include <Eigen/Eigenvalues> 

namespace fk {

template <class lattice_t>
struct configuration {
    const lattice_t& lattice;
    int_array_t f_config;
    double U, mu_c, mu_f;

    mutable real_matrix_t hamilt;
    mutable real_array_t cached_spectrum;
    mutable real_array_t cached_weights;
    mutable real_matrix_t cached_evecs;

    configuration(
        const lattice_t &lattice, double U, double mu_c, double mu_f):
        lattice(lattice),f_config(lattice.get_hopping_matrix().shape()[0]),
        U(U),mu_c(mu_c),mu_f(mu_f)
            { f_config()=0; }

    configuration(const configuration&) = default;
    configuration& operator= (const configuration& rhs) {U = rhs.U, mu_c = rhs.mu_c; mu_f = rhs.mu_f; 
                                                         f_config = rhs.f_config; 
                                                         hamilt = rhs.hamilt;
                                                         cached_weights = rhs.cached_weights;
                                                         cached_spectrum = rhs.cached_spectrum; 
                                                         cached_evecs = rhs.cached_evecs;
                                                         return (*this);};
    size_t get_nf() const;
    void randomize_f(triqs::mc_tools::random_generator &rnd, size_t nf = 0);

    real_matrix_t calc_hamiltonian();
    real_array_t  calc_spectrum();
    real_array_t  calc_spectrum_iter(real_array_t evals_guess, real_matrix_t evecs_guess, double beta);
};

template <class lattice_t>
inline size_t configuration<lattice_t>::get_nf() const
{
    return std::accumulate(f_config.begin(), f_config.end(), 0);
}

template <class lattice_t>
void configuration<lattice_t>::randomize_f(triqs::mc_tools::random_generator &rnd, size_t nf){
    if (!nf) nf = rnd(lattice.m_size);
    f_config()=0;
    for (size_t i=0; i<nf; ++i) {  
    size_t ind = rnd(lattice.m_size);
    while (f_config(ind)==1) ind = rnd(lattice.m_size);
    f_config(ind) = 1; 
    };
}

template <class lattice_t>
inline real_matrix_t configuration<lattice_t>::calc_hamiltonian()
{
    hamilt = lattice.get_hopping_matrix();
    for (size_t i=0; i<lattice.m_size; ++i) hamilt(i,i)+= -mu_c + U*f_config(i); // unoptimized
    return hamilt;
}

template <class lattice_t>
inline real_array_t configuration<lattice_t>::calc_spectrum()
{
    real_matrix_t h(hamilt);
    real_matrix_view_t h2(h);
    size_t size = hamilt.shape()[0]; 
    real_array_t evals(size);
    //std::tie(cached_spectrum, cached_evecs) = triqs::arrays::linalg::eigenelements(hamilt(),true);
    cached_spectrum = triqs::arrays::linalg::eigenvalues(h2);
    return cached_spectrum;
}

template <class lattice_t>
inline real_array_t configuration<lattice_t>::calc_spectrum_iter(real_array_t evals_guess, real_matrix_t evecs_guess, double beta)
{
/*
    Eigen::Map<EMatrixType<double>> map1 (&h.storage()[0], size, size);
    Eigen::SelfAdjointEigenSolver<EMatrixType<double>> s(map1);
    auto evals_eigen = s.eigenvalues();
    std::copy(evals_eigen.data(),evals_eigen.data()+size,evals.begin());
*/

}




} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_CONFIGURATION_HPP_

