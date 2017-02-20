#include "fk_mc/lattice/hypercubic.hpp"

namespace fk { 



template <size_t D>
hypercubic_lattice<D>::hypercubic_lattice(size_t lattice_size):
    abstract_lattice(sparse_m(boost::math::pow<D>(lattice_size), boost::math::pow<D>(lattice_size))),
    ft_pi_array_(msize())
{
    hopping_m_.reserve(Eigen::ArrayXi::Constant(msize(),D*2));
    dims_.fill(lattice_size);
    for (size_t i=0; i< msize(); i++) {
        auto pos = index_to_pos(i); 
        int v = 1;
        for (int p : pos) v*=((p%2)*2-1);
        ft_pi_array_[i]=v;
        }; 
};

template <size_t D>
std::array<int, D> hypercubic_lattice<D>::index_to_pos(size_t index) const
{
    std::array<int, D> out;
    for (int i=D-1; i>=0; i--) {
        out[i]=index%dims_[i];
        index/=dims_[i];
    };
    return out;
}

template <size_t D>
size_t hypercubic_lattice<D>::pos_to_index(std::array<int, D> pos) const
{
    size_t out=0;
    size_t mult = 1;
    for (int i=D-1; i>=0; i--) {
        out+=pos[i]*mult;
        mult*=dims_[i];
    };
    return out;
}

template <size_t D>
std::array<std::array<int, D>, 2*D> hypercubic_lattice<D>::neighbor_pos(std::array<int, D> pos) const {
    std::array<std::array<int, D>, 2*D> result;

    for (int i = 0; i < D; ++i) {
        result[2*i].fill(0);
        result[2*i+1].fill(0);
        result[2*i][i] = pos[i] == 0 ? dims_[i] - 1 : pos[i] - 1;
        result[2*i+1][i] = pos[i] == dims_[i] - 1 ? 0 : pos[i] + 1;
    }

    return result;
}

template <size_t D>
std::vector<size_t> hypercubic_lattice<D>::neighbor_index(size_t index) const {
    std::vector<size_t> result(2*D);
    auto positions = neighbor_pos(index_to_pos(index));
    std::transform(std::begin(positions), std::end(positions),
                   std::begin(result),
                   [this](const std::array<int, D> &p) { return pos_to_index(p); });
    return result;
}

template <size_t D>
int hypercubic_lattice<D>::FFT_pi(const Eigen::ArrayXi& in) const
{
    return (in*ft_pi_array_).sum();
}


template <size_t D>
typename hypercubic_lattice<D>::BZPoint hypercubic_lattice<D>::get_bzpoint(std::array<double, D> in) const
{
    return BZPoint(in, *this);
}

template <size_t D>
typename hypercubic_lattice<D>::BZPoint hypercubic_lattice<D>::get_bzpoint(size_t in) const
{
    return BZPoint(in, *this);
}
    
template <size_t D>
std::vector<typename hypercubic_lattice<D>::BZPoint> hypercubic_lattice<D>::get_all_bzpoints() const
{
    int npts=1; for (auto x:dims_) npts*=x;
    std::vector<BZPoint> out;
    out.reserve(npts);
    for (int i=0; i<npts; i++) out.push_back(this->get_bzpoint(i));
    return out;
}


template struct hypercubic_lattice<1>;
template struct hypercubic_lattice<2>;
template struct hypercubic_lattice<3>;
//template struct hypercubic_lattice<4>;


template <size_t D>
hypercubic_lattice<D>& fill_nearest_neighbors(hypercubic_lattice<D> &l, double t)
{
    for (size_t i=0; i< l.volume(); ++i) {
        for (int o=0; o< l.norbs(); ++o) {
            auto current_pos = l.index_to_pos(i);
            for (size_t n=0; n<D; ++n) {
                auto pos_l(current_pos), pos_r(current_pos);
                pos_l[n]=(current_pos[n]>0?current_pos[n]-1:l.dims()[n]-1);
                pos_r[n]=(current_pos[n]<l.dims()[n]-1?current_pos[n]+1:0);
                l.add_hopping(i, l.pos_to_index(pos_l), -1.0*t);
                l.add_hopping(i, l.pos_to_index(pos_r), -1.0*t);
            };
        };
    };
    return l;
}

template hypercubic_lattice<1>& fill_nearest_neighbors (hypercubic_lattice<1> &, double);
template hypercubic_lattice<2>& fill_nearest_neighbors (hypercubic_lattice<2> &, double);
template hypercubic_lattice<3>& fill_nearest_neighbors (hypercubic_lattice<3> &, double);

hypercubic_lattice<2>& fill_triangular(hypercubic_lattice<2> &l, double t, double tp)
{
    fill_nearest_neighbors(l, t);
    for (size_t i=0; i<l.volume(); ++i) {
        for (int o=0; o< l.norbs(); ++o) {
            auto current_pos = l.index_to_pos(i);
            auto pos_l(current_pos), pos_r(current_pos);
            for (size_t n=0; n<2; ++n) {
                pos_l[n]=(current_pos[n]>0?current_pos[n]-1:l.dims()[n]-1);
                pos_r[n]=(current_pos[n]<l.dims()[n]-1?current_pos[n]+1:0);
                };

            l.add_hopping(i, l.pos_to_index(pos_l), -1.0*tp);
            l.add_hopping(i, l.pos_to_index(pos_r), -1.0*tp);
            };
        }

    return l;
}




hypercubic_lattice<2>& fill_honeycomb(hypercubic_lattice<2> &l, double t)
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
        auto dims = l.dims();

        if (l.dims()[X] != l.dims()[Y]) throw std::logic_error("X dim != Y dim");
        if (l.dims()[X] %2 != 0 || l.dims()[Y] %2 != 0) throw std::logic_error("Need even size");

        bool subA = true;
        for (size_t i=0; i<l.volume(); ++i) {
            auto current_pos = l.index_to_pos(i);
            auto pos_l(current_pos), pos_r(current_pos),pos_u(current_pos),pos_d(current_pos);

            pos_l[X]=(current_pos[X]>0?current_pos[X]-1:dims[X]-1);
            pos_r[X]=(current_pos[X]<dims[X]-1?current_pos[X]+1:0);
            pos_d[Y]=(current_pos[Y]>0?current_pos[Y]-1:dims[Y]-1);
            pos_u[Y]=(current_pos[Y]<dims[Y]-1?current_pos[Y]+1:0);

            l.add_hopping(i,l.pos_to_index(pos_l), -1.0*t);
            l.add_hopping(i,l.pos_to_index(pos_r), -1.0*t);
            if (subA)
                l.add_hopping(i,l.pos_to_index(pos_u), -1.0*t);
            else
                l.add_hopping(i,l.pos_to_index(pos_d), -1.0*t);
            subA = !subA;
        };
   return l;
}


} // end of namespace fk
