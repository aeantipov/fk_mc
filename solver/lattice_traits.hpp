#ifndef __FK_LATTICE_TRAITS_HPP
#define __FK_LATTICE_TRAITS_HPP

#include "common.hpp"
#include <Eigen/SparseCore>

using namespace triqs;

namespace fk {

/** A class that represents a lattice and generates the hopping matrix 
 * of a hypercubic lattice in different dimensions. */ 
template <size_t D>
struct square_lattice_traits 
{
    typedef triqs::arrays::matrix<double>  matrix_t;
    typedef triqs::arrays::matrix_view<double>  matrix_view_t;
    typedef Eigen::SparseMatrix<double> sparse_m;

    static constexpr size_t Ndim = D;
    std::array<size_t, D> dims;
    size_t m_size;

    matrix_t hopping_m;
    sparse_m hopping_sparse;
    size_t nonzero_elems = 0;

    std::array<size_t, D> index_to_pos(size_t index);
    size_t pos_to_index(std::array<size_t, D> pos);

    //triqs::arrays::array_view<double,D> matrix_view ( real_array_view_t in );
    //real_array_view_t flatten(triqs::arrays::array_view<double,D> in);

    square_lattice_traits(size_t lattice_size):
        m_size(boost::math::pow<D>(lattice_size)),
        hopping_m(m_size,m_size)
        ,hopping_sparse(m_size,m_size)
    {
            hopping_sparse.reserve(Eigen::VectorXi::Constant(m_size,D*2));
            dims.fill(lattice_size);
            hopping_m()=0;
    };

    void fill(double t);

    matrix_view_t get_hopping_matrix() const {return hopping_m;};
};

/** Specification of square_lattice_traits for the triangular lattice. */
struct triangular_lattice_traits : square_lattice_traits<2> 
{
    using square_lattice_traits<2>::dims;
    using square_lattice_traits<2>::m_size;
    using square_lattice_traits<2>::hopping_m;
    triangular_lattice_traits(size_t lattice_size):square_lattice_traits<2>(lattice_size){
        hopping_sparse.reserve(Eigen::VectorXi::Constant(m_size,6));
        };
    void fill(double t, double t_p);
};

}; // end of namespace FK

#endif // #ifndef __FK_LATTICE_TRAITS_HPP
