#include "triangular.hpp"

namespace fk { 

void triangular_lattice::fill(double t, double t_p)
{
    hypercubic_lattice<2>::fill(t);
    for (size_t i=0; i<m_size_; ++i) {
        auto current_pos = index_to_pos(i);
        auto pos_l(current_pos), pos_r(current_pos);
        for (size_t n=0; n<2; ++n) {
            pos_l[n]=(current_pos[n]>0?current_pos[n]-1:dims[n]-1);
            pos_r[n]=(current_pos[n]<dims[n]-1?current_pos[n]+1:0);
            };

        hopping_m.insert(i,pos_to_index(pos_l)) = -1.0*t_p;
        hopping_m.insert(i,pos_to_index(pos_r)) = -1.0*t_p;
        };
}

} // end of namespace fk
