#ifndef __FK_LATTICE_TRAITS_HPP
#define __FK_LATTICE_TRAITS_HPP

#include "common.hpp"
#include <Eigen/SparseCore>
#include <fftw3.h>

using namespace triqs;

namespace fk {

struct lattice_base {
    typedef Eigen::MatrixXd dense_m;
    typedef Eigen::SparseMatrix<double> sparse_m;
    int get_msize() const {return m_size_;}; 
    sparse_m hopping_m;
    lattice_base(sparse_m &&in); // hopping_m(in),m_size_(hopping_m.rows())
protected:
    size_t m_size_;
};


/** A class that represents a lattice and generates the hopping matrix 
 * of a hypercubic lattice in different dimensions. */ 
template <size_t D>
struct hypercubic_lattice : lattice_base
{

    std::array<size_t, D> index_to_pos(size_t index);
    size_t pos_to_index(std::array<size_t, D> pos);

    static constexpr size_t Ndim = D;
    std::array<int, D> dims;
    using lattice_base::m_size_;
    using lattice_base::hopping_m;

    Eigen::VectorXcd FFT(Eigen::VectorXcd in, int direction);
    //triqs::arrays::array_view<double,D> matrix_view ( real_array_view_t in );
    //real_array_view_t flatten(triqs::arrays::array_view<double,D> in);

    hypercubic_lattice(size_t lattice_size);
        
    void fill(double t);
};

/** Specification of hypercubic_lattice for the triangular lattice. */
struct triangular_lattice : hypercubic_lattice<2> 
{
    using hypercubic_lattice<2>::dims;
    using hypercubic_lattice<2>::m_size_;
    using hypercubic_lattice<2>::hopping_m;
    triangular_lattice(size_t lattice_size):hypercubic_lattice<2>(lattice_size){
        hopping_m.reserve(Eigen::VectorXi::Constant(m_size_,6));
        };
    void fill(double t, double t_p);
};

}; // end of namespace FK

#endif // #ifndef __FK_LATTICE_TRAITS_HPP
