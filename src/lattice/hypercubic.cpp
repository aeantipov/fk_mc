#include "fk_mc/lattice/hypercubic.hpp"

namespace fk { 

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

template <size_t D>
hypercubic_lattice<D>::hypercubic_lattice(size_t lattice_size):
    lattice_base(sparse_m(boost::math::pow<D>(lattice_size), boost::math::pow<D>(lattice_size))),
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

template hypercubic_lattice<1>& fill_nearest_neighbors (hypercubic_lattice<1> &, double);
template hypercubic_lattice<2>& fill_nearest_neighbors (hypercubic_lattice<2> &, double);
template hypercubic_lattice<3>& fill_nearest_neighbors (hypercubic_lattice<3> &, double);

} // end of namespace fk
