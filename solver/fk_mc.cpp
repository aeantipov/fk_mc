#include "fk_mc.hxx"

#include "lattice/hypercubic.hpp"
#include "lattice/triangular.hpp"
#include "lattice/chain.hpp"
#include "lattice/honeycomb.hpp"

namespace fk {

template class fk_mc<triangular_lattice>; 
template class fk_mc<hypercubic_lattice<1>>; 
template class fk_mc<hypercubic_lattice<2>>; 
template class fk_mc<chain_lattice>; 
template class fk_mc<honeycomb_lattice>; 

void observables_t::reserve(int n) 
{ 
    energies.reserve(n);
    d2energies.reserve(n);
    nf0.reserve(n);
    nfpi.reserve(n);
    spectrum.reserve(n); 
    spectrum_history.reserve(n); // L^D x n_measures size
    ipr_history.reserve(n); // L^D x n_measures size
    focc_history.reserve(n);     // L^D x n_measures size
    nq_history.reserve(n);       // nqpts x n_measures size
    fsuscq_history.reserve(n);   // nqpts x n_measures size
    eigenfunctions_history.reserve(n);   // L^D * L^D * n_measures size (huge)
}

template<typename T>
void auto_merge(T& in, T& out)
{
    if (out.size() == 0 && in.size() == 0 ) { return; }
    if (in.size() == 0) { std::swap(in, out); return; } 
    if (out.size() > 0) std::move(out.begin(), out.end(), std::back_inserter(in));
}

template<typename T>
void auto_merge_vv(T& in, T& out)
{
    if (out.size() == 0 && in.size() == 0 ) { return; }
    if (in.size() == 0) { std::swap(in, out); return; } 
    if (out.size() > 0) { for (int i=0; i<out.size(); i++) auto_merge(in[i], out[i]); } 
}

void observables_t::merge(observables_t& rhs)
{
    if (rhs.nfpi.size()==0) return;
    if (this->nfpi.size()==0) { std::swap(*this, rhs); return; }
    auto_merge(energies, rhs.energies);
    auto_merge(d2energies, rhs.d2energies);
    auto_merge(nf0, rhs.nf0);
    auto_merge(nfpi, rhs.nfpi);
    auto_merge(spectrum, rhs.spectrum);
    auto_merge(eigenfunctions_history, rhs.eigenfunctions_history);

    auto_merge_vv(spectrum_history, rhs.spectrum_history);
    auto_merge_vv(ipr_history, rhs.ipr_history);
    auto_merge_vv(focc_history, rhs.focc_history);
    auto_merge_vv(nq_history, rhs.nq_history);
    auto_merge_vv(fsuscq_history, rhs.fsuscq_history);
} 


} // end of namespace FK
