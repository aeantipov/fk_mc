#include "lattice.hpp"


namespace fk { 

lattice_base::lattice_base(sparse_m &&in):
   hopping_m(in),
   m_size_(hopping_m.rows()) 
{
    if (hopping_m.rows() != hopping_m.cols() || hopping_m.rows() == 0) TRIQS_RUNTIME_ERROR << "Failed to initalize lattice. ";
}

template <size_t D>
hypercubic_lattice<D>::hypercubic_lattice(size_t lattice_size):
    lattice_base(sparse_m(boost::math::pow<D>(lattice_size), boost::math::pow<D>(lattice_size))),
    ft_pi_array_(m_size_)
{
    hopping_m.reserve(Eigen::ArrayXi::Constant(m_size_,D*2));
    dims.fill(lattice_size);
    for (size_t i=0; i<m_size_; i++) { 
        auto pos = index_to_pos(i); 
        int v = 1;
        for (int p : pos) v*=((p%2)*2-1);
        ft_pi_array_[i]=v;
        }; 
};

template <size_t D>
std::array<size_t, D> hypercubic_lattice<D>::index_to_pos(size_t index) const
{
    std::array<size_t, D> out;
    for (int i=D-1; i>=0; i--) {
        out[i]=index%dims[i];
        index/=dims[i];
    };
    return out;
}

template <size_t D>
size_t hypercubic_lattice<D>::pos_to_index(std::array<size_t, D> pos) const
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
Eigen::ArrayXcd hypercubic_lattice<D>::FFT(Eigen::ArrayXcd in, int direction) const
{
    Eigen::ArrayXcd out(in);

    fftw_plan p;
    p = fftw_plan_dft(D, dims.data(),
                         reinterpret_cast<fftw_complex*>( in.data()), 
                         reinterpret_cast<fftw_complex*>( out.data()), 
                         direction, FFTW_ESTIMATE); 
    fftw_execute(p);

    double norm=1.0;
    for (auto x:dims) norm*=x;
    if (direction == FFTW_BACKWARD) out/=norm;
    return out;

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
    int npts=1; for (auto x:dims) npts*=x;
    std::vector<BZPoint> out;
    out.reserve(npts);
    for (int i=0; i<npts; i++) out.push_back(this->get_bzpoint(i));
    return out;
}

template <size_t D>
void hypercubic_lattice<D>::fill(double t)
{
    for (size_t i=0; i<m_size_; ++i) {
        auto current_pos = index_to_pos(i);
        for (size_t n=0; n<D; ++n) {
            auto pos_l(current_pos), pos_r(current_pos);
            pos_l[n]=(current_pos[n]>0?current_pos[n]-1:dims[n]-1);
            pos_r[n]=(current_pos[n]<dims[n]-1?current_pos[n]+1:0);
            hopping_m.insert(i,pos_to_index(pos_l)) = -1.0*t;
            hopping_m.insert(i,pos_to_index(pos_r)) = -1.0*t;
        }; 
    };
}


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

template struct hypercubic_lattice<1>;
template struct hypercubic_lattice<2>;
template struct hypercubic_lattice<3>;
//template struct hypercubic_lattice<4>;

} // end of namespace fk
