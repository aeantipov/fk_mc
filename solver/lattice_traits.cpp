#include "lattice_traits.hpp"

namespace fk { 

template <size_t D>
inline std::array<size_t, D> square_lattice_traits<D>::index_to_pos(size_t index)
{
    std::array<size_t, D> out;
    for (int i=D-1; i>=0; i--) {
        out[i]=index%dims[i];
        index/=dims[i];
    };
    return out;
}

template <size_t D>
inline size_t square_lattice_traits<D>::pos_to_index(std::array<size_t, D> pos)
{
    size_t out=0;
    size_t mult = 1;
    for (int i=D-1; i>=0; i--) {
        out+=pos[i]*mult;
        mult*=dims[i];
    };
    return out;
}

template <size_t D>
void square_lattice_traits<D>::fill(double t)
{
    for (size_t i=0; i<m_size; ++i) {
        auto current_pos = index_to_pos(i);
        for (size_t n=0; n<D; ++n) {
            auto pos_l(current_pos), pos_r(current_pos);
            pos_l[n]=(current_pos[n]>0?current_pos[n]-1:dims[n]-1);
            pos_r[n]=(current_pos[n]<dims[n]-1?current_pos[n]+1:0);
            hopping_m(i,pos_to_index(pos_l)) = t;
            hopping_m(i,pos_to_index(pos_r)) = t;
        }; 
    };
}


void triangular_lattice_traits::fill(double t, double t_p)
{
    square_lattice_traits<2>::fill(t);
    for (size_t i=0; i<m_size; ++i) {
        auto current_pos = index_to_pos(i);
        auto pos_l(current_pos), pos_r(current_pos);
        for (size_t n=0; n<2; ++n) {
            pos_l[n]=(current_pos[n]>0?current_pos[n]-1:dims[n]-1);
            pos_r[n]=(current_pos[n]<dims[n]-1?current_pos[n]+1:0);
            };

        hopping_m(i,pos_to_index(pos_l)) = t_p;
        hopping_m(i,pos_to_index(pos_r)) = t_p;

        };
}

template struct square_lattice_traits<1>;
template struct square_lattice_traits<2>;
template struct square_lattice_traits<3>;
template struct square_lattice_traits<4>;

} // end of namespace fk
