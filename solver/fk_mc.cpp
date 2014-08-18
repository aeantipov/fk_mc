#include "fk_mc.hxx"

#include "lattice/hypercubic.hpp"
#include "lattice/triangular.hpp"
#include "lattice/chain.hpp"

namespace fk {

template class fk_mc<triangular_lattice>; 
template class fk_mc<hypercubic_lattice<1>>; 
template class fk_mc<hypercubic_lattice<2>>; 
template class fk_mc<chain_lattice>; 

void observables_t::reserve(int n) 
{ 
    energies.reserve(n);
    d2energies.reserve(n);
    nf0.reserve(n);
    nfpi.reserve(n);
    spectrum_history.reserve(n); // L^D x n_measures size
    ipr_history.reserve(n); // L^D x n_measures size
    focc_history.reserve(n);     // L^D x n_measures size
    nq_history.reserve(n);       // nqpts x n_measures size
    fsuscq_history.reserve(n);   // nqpts x n_measures size
}



} // end of namespace FK
