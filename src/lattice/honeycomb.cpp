#include "fk_mc/lattice/honeycomb.hpp"

namespace fk { 

void honeycomb_lattice::fill(double t)
{
    /** The rule for filling the lattice is the following:
 y ^
   | -B-A-B-A-
   | -|---|---
   | -A-B-A-B-
   | ---|---|-
   | -B-A-B-A-
   | -|---|---
   | -A-B-A-B-
   |___________> x
    A always connects to B, and connection always 
    go up from A sites and down from B
    x will be the second axis, y - the first
    */

    constexpr static int X = 1;
    constexpr static int Y = 0;

    if (dims[X] != dims[Y]) throw std::logic_error("X dim != Y dim"); 
    if (dims[X] %2 != 0 || dims[Y] %2 != 0) throw std::logic_error("Need even size");

    bool subA = true;
    for (size_t i=0; i<m_size_; ++i) {
        auto current_pos = index_to_pos(i);
        auto pos_l(current_pos), pos_r(current_pos),pos_u(current_pos),pos_d(current_pos);

        pos_l[X]=(current_pos[X]>0?current_pos[X]-1:dims[X]-1);
        pos_r[X]=(current_pos[X]<dims[X]-1?current_pos[X]+1:0);
        pos_d[Y]=(current_pos[Y]>0?current_pos[Y]-1:dims[Y]-1);
        pos_u[Y]=(current_pos[Y]<dims[Y]-1?current_pos[Y]+1:0);

        hopping_m_.insert(i,pos_to_index(pos_l)) = -1.0*t;
        hopping_m_.insert(i,pos_to_index(pos_r)) = -1.0*t;
        if (subA)  
            hopping_m_.insert(i,pos_to_index(pos_u)) = -1.0*t;
        else 
            hopping_m_.insert(i,pos_to_index(pos_d)) = -1.0*t;
        subA = !subA; 
    };
}

} // end of namespace fk
