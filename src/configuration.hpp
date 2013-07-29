#ifndef __FK_MC_CONFIGURATION_HPP_
#define __FK_MC_CONFIGURATION_HPP_

#include "common.hpp"
#include <triqs/mc_tools/random_generator.hpp>
#include <numeric>
#include <triqs/arrays/linalg/eigenelements.hpp>

namespace fk {

template <class lattice_t>
struct configuration {
    const lattice_t& lattice;
    int_array_t f_config;
    double U, mu_c, mu_f;
    
    mutable real_array_t cached_spectrum;

    configuration(
        const lattice_t &lattice,
        double U,
        double mu_c,
        double mu_f
                ):
        lattice(lattice),f_config(lattice.get_hopping_matrix().shape()[0]),
        U(U),mu_c(mu_c),mu_f(mu_f)
        {
        f_config()=0;
        }

    configuration(const configuration&) = default;
    configuration& operator= (const configuration& rhs) {U = rhs.U, mu_c = rhs.mu_c; mu_f = rhs.mu_f; 
                                                         f_config = rhs.f_config; cached_spectrum = rhs.cached_spectrum; return (*this);};
    size_t get_nf() const;
    void randomize_f(triqs::mc_tools::random_generator &rnd, size_t nf = 0);

    real_matrix_t get_hamiltonian() const;
    real_array_t get_spectrum() const;
};

template <class lattice_t>
inline size_t configuration<lattice_t>::get_nf() const
{
    return std::accumulate(f_config.begin(), f_config.end(), 0);
}

template <class lattice_t>
void configuration<lattice_t>::randomize_f(triqs::mc_tools::random_generator &rnd, size_t nf){
    if (!nf) nf = get_nf();
    for (size_t i=0; i<nf; ++i) {  
    size_t ind = rnd(lattice.m_size);
    while (f_config(ind)==1) ind = rnd(lattice.m_size);
    f_config(ind) = 1; 
    };
}

template <class lattice_t>
real_matrix_t configuration<lattice_t>::get_hamiltonian() const
{
    real_matrix_t T1(lattice.get_hopping_matrix());
    for (size_t i=0; i<lattice.m_size; ++i) T1(i,i)+= -mu_c + U*f_config(i); // unoptimized
    return T1;
}

template <class lattice_t>
real_array_t configuration<lattice_t>::get_spectrum() const
{
    auto evals = triqs::arrays::linalg::eigenvalues(real_matrix_view_t(get_hamiltonian())); 
    cached_spectrum = evals;
    return evals;
}

// Free functions
inline real_array_t density_matrix_c(double beta, real_array_t spectrum, double offset_energy)
{
    double exp_offset = exp(beta*offset_energy);
    auto F = triqs::arrays::map(std::function<double(double)>( [beta,offset_energy,exp_offset](double E){return 1.0+exp(-beta*(E-offset_energy))/exp_offset;} ));
    return F(spectrum);
}



} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_CONFIGURATION_HPP_

