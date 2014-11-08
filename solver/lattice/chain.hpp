#pragma once

#include "lattice/hypercubic.hpp"

namespace fk { 

/** Specification of hypercubic_lattice for the triangular lattice. */
struct chain_lattice : hypercubic_lattice<1> 
{
    using hypercubic_lattice<1>::Ndim;
    using hypercubic_lattice<1>::dims;
    using hypercubic_lattice<1>::m_size_;
    using hypercubic_lattice<1>::hopping_m_;
    chain_lattice(size_t lattice_size):hypercubic_lattice<1>(lattice_size){
        hopping_m_.reserve(Eigen::ArrayXi::Constant(m_size_,2));
        };
    void fill(double t, double eta, double delta);
};

} // end of namespace fk
