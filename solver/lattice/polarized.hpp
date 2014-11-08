#pragma once

#include "lattice/hypercubic.hpp"

namespace fk { 

/** Polarized hypercubic lattice. */
template <size_t D>
struct polarized_lattice : hypercubic_lattice<D> 
{
    typedef hypercubic_lattice<D> base;
    using base::Ndim;
    using base::dims;
    using base::m_size_;
    using base::hopping_m_;
    polarized_lattice(size_t lattice_size):base(lattice_size){}
    void fill(double t1, double t2);
    void fill(double t1) { this->fill(t1,t1); }
};

template <size_t D>
inline void polarized_lattice<D>::fill(double t1, double t2)
{
    for (size_t i=0; i<m_size_; ++i) {
        auto pos = this->index_to_pos(i);
        int is_even = 0;
        for (size_t n=0; n<Ndim; ++n) is_even += pos[n]%2;
        is_even=1-is_even%2;
        FKDEBUG(pos << " even? : " << is_even);
        auto pos_l(pos), pos_r(pos);
        for (size_t n=0; n<Ndim; ++n) {
            pos_l[n]=(pos[n]>0?pos[n]-1:dims[n]-1);
            pos_r[n]=(pos[n]<dims[n]-1?pos[n]+1:0);
            if (n==0) {
                // introduce polarization in x axis
                if (is_even) { 
                    hopping_m_.insert(i,this->pos_to_index(pos_l)) = -1.0*t1;
                    hopping_m_.insert(i,this->pos_to_index(pos_r)) = -1.0*t2;
                }
                else {
                    hopping_m_.insert(i,this->pos_to_index(pos_l)) = -1.0*t2;
                    hopping_m_.insert(i,this->pos_to_index(pos_r)) = -1.0*t1;
                    }
                }
            else { 
                // regular hopping in y direction 
                hopping_m_.insert(i,this->pos_to_index(pos_l)) = -1.0*t1;
                hopping_m_.insert(i,this->pos_to_index(pos_r)) = -1.0*t1;
            }
           // FKDEBUG(pos << " -> " << pos_l << " : " << hopping_m_.coeffRef(i,this->pos_to_index(pos_l)));
           // FKDEBUG(pos << " -> " << pos_r << " : " << hopping_m_.coeffRef(i,this->pos_to_index(pos_r)));

        };
        };
}



} // end of namespace fk
