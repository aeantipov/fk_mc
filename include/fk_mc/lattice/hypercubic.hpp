#pragma once 

#include <Eigen/SparseCore>
#include <fftw3.h>

#include "common.hpp"
#include "lattice.hpp"

//using namespace triqs;

namespace fk {

/** A class that represents a lattice and generates the hopping matrix 
 * of a hypercubic lattice in different dimensions. */ 
template <size_t D>
struct hypercubic_lattice : lattice_base
{

    std::array<int, D> index_to_pos(size_t index) const;
    size_t pos_to_index(std::array<int, D> pos) const;

    struct BZPoint { std::array<double, D> val_; 
                     size_t ind_; 
                     explicit operator int() const {return ind_;};
                     explicit operator size_t() const {return ind_;};
                     operator std::array<double, D> () const {return val_;};
                   protected:
                     BZPoint(size_t ind, const hypercubic_lattice<D> &l):ind_(ind){ 
                        auto pos = l.index_to_pos(ind_); 
                        for (size_t i=0; i<D; i++) { val_[i] = double(pos[i])*2.0*PI/l.dims[i]; };
                        };
                     BZPoint(std::array<double, D> v, const hypercubic_lattice<D> &l, double prec = 1e-12){ 
                        std::array<int, D> pos;
                        for (size_t i=0; i<D; i++) { pos[i] = std::round(v[i]/2.0/PI*l.dims[i]); };
                        size_t i = l.pos_to_index(pos); 
                        (*this) = BZPoint(i,l); 
                        double diff = 0; for (size_t j=0; j<D; j++) diff+=std::abs(v[j]-val_[j]);
                        if (diff>=prec) TRIQS_RUNTIME_ERROR<<"Couldn't find right BZPoint.";
                        };
                    friend class hypercubic_lattice<D>;

                     friend std::ostream& operator<< (std::ostream& out, BZPoint q){ out << "{"; for (double x:q.val_) out << x << " "; out << "}"; return out;};
                   };

    BZPoint get_bzpoint(std::array<double, D> in) const;
    BZPoint get_bzpoint(size_t in) const;
    std::vector<BZPoint> get_all_bzpoints() const;

    static constexpr size_t Ndim = D;
    std::array<int, D> dims;
    using lattice_base::m_size_;
    using lattice_base::hopping_m_;

    template <typename M> M FFT(M in, int direction) const;
    int FFT_pi(const Eigen::ArrayXi& in) const;
    //triqs::arrays::array_view<double,D> matrix_view ( real_array_view_t in );
    //real_array_view_t flatten(triqs::arrays::array_view<double,D> in);

    hypercubic_lattice(size_t lattice_size);
    hypercubic_lattice(hypercubic_lattice const& rhs):lattice_base(rhs), dims(rhs.dims), ft_pi_array_(rhs.ft_pi_array_){} 
        
    void fill(double t);
//    protected:
        Eigen::ArrayXi ft_pi_array_;
};

template <size_t D>
template <typename M> 
M hypercubic_lattice<D>::FFT(M in, int direction) const
{
    M out(in);
    out.setZero();

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

}; // end of namespace FK

