#include "configuration.hpp"

namespace fk {

size_t configuration_t::get_nf() const
{
    return std::accumulate(f_config.data(), f_config.data()+lattice.get_msize(), 0);
}

void configuration_t::randomize_f(triqs::mc_tools::random_generator &rnd, size_t nf){
    if (!nf) nf = rnd(lattice.get_msize());
    f_config.setZero();
    for (size_t i=0; i<nf; ++i) {  
    size_t ind = rnd(lattice.get_msize());
    while (f_config(ind)==1) ind = rnd(lattice.get_msize());
    f_config(ind) = 1; 
    };
}

typename configuration_t::sparse_m configuration_t::calc_hamiltonian()
{
    hamilt.reserve(lattice.hopping_m.nonZeros() + lattice.get_msize());
    MY_DEBUG("!");
    hamilt = lattice.hopping_m;
    MY_DEBUG(lattice.hopping_m << " " << lattice.get_msize());
    for (size_t i=0; i<lattice.get_msize(); ++i) hamilt.coeffRef(i,i)+= -mu_c + U*f_config(i); // unoptimized
    MY_DEBUG("!!");
    return hamilt;
}


typename configuration_t::real_array_t configuration_t::calc_spectrum()
{
    dense_m h(hamilt);
    Eigen::SelfAdjointEigenSolver<dense_m> s(h,Eigen::EigenvaluesOnly);
    cached_spectrum = s.eigenvalues();
    std::sort (cached_spectrum.data(), cached_spectrum.data()+cached_spectrum.size());  

    cached_weights.resize(cached_spectrum.size());

    double e0 = cached_spectrum[0];
    double weight0 = exp(-beta*e0);
    for (size_t i=0; i<cached_spectrum.size(); ++i) { 
        cached_weights(i) = exp(-beta*(cached_spectrum(i)-e0)); 
        };

    return cached_spectrum;
}

/* Obsolete, clean 
typename configuration_t::real_array_t configuration_t::calc_spectrum_arpack()
{

    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<sparse_m> s(hamilt,n_calc_evals,"SA",Eigen::EigenvaluesOnly);
    cached_spectrum = s.eigenvalues().reverse();
    return cached_spectrum;
}
*/





} // end of namespace fk
