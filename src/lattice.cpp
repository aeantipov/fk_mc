#include "fk_mc/lattice.hpp"

namespace fk { 

template <typename T> inline T myconj(T in);
template <> inline std::complex<double> myconj<std::complex<melem_type>> (std::complex<melem_type> in) { return std::conj(in); }
template <> inline double myconj<double> (double in) { return in; }

abstract_lattice& abstract_lattice::add_hopping(size_t l, size_t r, melem_type v, bool symmetrize)
{ 
    hopping_m_.insert(l,r) = v;
    if (symmetrize) { 
        hopping_m_.insert(r,l) = myconj(v); 
        }
    return (*this);
}

} // end of namespace fk
