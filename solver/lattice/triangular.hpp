#pragma once

#include "lattice/hypercubic.hpp"

namespace fk { 

/** Specification of hypercubic_lattice for the triangular lattice. */
struct triangular_lattice : hypercubic_lattice<2> 
{
    using hypercubic_lattice<2>::Ndim;
    using hypercubic_lattice<2>::dims;
    using hypercubic_lattice<2>::m_size_;
    using hypercubic_lattice<2>::hopping_m;
    triangular_lattice(size_t lattice_size):hypercubic_lattice<2>(lattice_size){
        hopping_m.reserve(Eigen::ArrayXi::Constant(m_size_,6));
        };
    void fill(double t, double t_p);
};

} // end of namespace fk
