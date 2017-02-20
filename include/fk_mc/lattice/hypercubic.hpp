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
class hypercubic_lattice : public abstract_lattice
{
    using abstract_lattice::hopping_m_;
public:
    static constexpr size_t Ndim = D;

    std::array<int, D> index_to_pos(size_t index) const;
    size_t pos_to_index(std::array<int, D> pos) const;
    /// The nearest neighbor positions to the given position.
    std::array<std::array<int, D>, 2*D> neighbor_pos(std::array<int, D> pos) const;
    /// The nearest neighbor indices to the given index.
    std::vector<size_t> neighbor_index(size_t index) const override;
    size_t ndim() const override { return Ndim; }
    struct BZPoint { std::array<double, D> val_; 
                     size_t ind_; 
                     explicit operator int() const {return ind_;};
                     explicit operator size_t() const {return ind_;};
                     operator std::array<double, D> () const {return val_;};
                   protected:
                     BZPoint(size_t ind, const hypercubic_lattice<D> &l):ind_(ind){ 
                        auto pos = l.index_to_pos(ind_); 
                        for (size_t i=0; i<D; i++) { val_[i] = double(pos[i])*2.0*PI/l.dims_[i]; };
                        };
                     BZPoint(std::array<double, D> v, const hypercubic_lattice<D> &l, double prec = 1e-12){ 
                        std::array<int, D> pos;
                        for (size_t i=0; i<D; i++) { pos[i] = std::round(v[i]/2.0/PI*l.dims_[i]); };
                        size_t i = l.pos_to_index(pos); 
                        (*this) = BZPoint(i,l); 
                        double diff = 0; for (size_t j=0; j<D; j++) diff+=std::abs(v[j]-val_[j]);
                        if (diff>=prec) FKMC_ERROR<<"Couldn't find right BZPoint.";
                        };
                    friend class hypercubic_lattice<D>;

                    friend std::ostream& operator<< (std::ostream& out, BZPoint q){ out << "{"; for (double x:q.val_) out << x << " "; out << "}"; return out;};
                   };

    BZPoint get_bzpoint(std::array<double, D> in) const;
    BZPoint get_bzpoint(size_t in) const;
    std::vector<BZPoint> get_all_bzpoints() const;
    template <typename M> M FFT(M in, int direction) const;
    int FFT_pi(const Eigen::ArrayXi& in) const;
    //triqs::arrays::array_view<double,D> matrix_view ( real_array_view_t in );
    //real_array_view_t flatten(triqs::arrays::array_view<double,D> in);

    hypercubic_lattice(size_t lattice_size);
    hypercubic_lattice(hypercubic_lattice const& rhs):abstract_lattice(rhs), dims_(rhs.dims_), ft_pi_array_(rhs.ft_pi_array_){}
        
    std::array<int, D> const& dims() const { return dims_; }

protected:
    std::array<int, D> dims_;
    Eigen::ArrayXi ft_pi_array_;
};

template <size_t D>
template <typename M> 
M hypercubic_lattice<D>::FFT(M in, int direction) const
{
    M out(in);
    out.setZero();

    fftw_plan p;
    p = fftw_plan_dft(D, dims_.data(),
                         reinterpret_cast<fftw_complex*>( in.data()), 
                         reinterpret_cast<fftw_complex*>( out.data()), 
                         direction, FFTW_ESTIMATE); 
    fftw_execute(p);

    double norm=1.0;
    for (auto x:dims_) norm*=x;
    if (direction == FFTW_BACKWARD) out/=norm;
    return out;
}

template <size_t D>
hypercubic_lattice<D>& fill_nearest_neighbors(hypercubic_lattice<D> &l, double t);

hypercubic_lattice<2>& fill_triangular(hypercubic_lattice<2> &l, double t, double tp);
hypercubic_lattice<2>& fill_honeycomb(hypercubic_lattice<2> &l, double t);

}; // end of namespace FK

