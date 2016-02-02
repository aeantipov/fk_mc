#pragma once

#include "lattice/hypercubic.hpp"

namespace fk { 

/** Specification of hypercubic_lattice for the honeycomb lattice. 
    The ``brickwall'' representation of the honeycomb lattice is used.
 */
struct honeycomb_lattice : hypercubic_lattice<2> 
{
    using hypercubic_lattice<2>::Ndim;
    using hypercubic_lattice<2>::dims;
    using hypercubic_lattice<2>::m_size_;
    using hypercubic_lattice<2>::hopping_m_;
    honeycomb_lattice(size_t lattice_size):hypercubic_lattice<2>(lattice_size){
        hopping_m_.reserve(Eigen::ArrayXi::Constant(m_size_,3));
        };
    void fill(double t);
};

} // end of namespace fk
