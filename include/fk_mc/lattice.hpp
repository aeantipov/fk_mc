#pragma once 

#include "common.hpp"
#include <Eigen/SparseCore>
#include <fftw3.h>

namespace fk {

/** A class to represent a lattice of a finite volume. 
 *  It defines and provides the tight-binding "hopping" matrix in real space.
 */
struct lattice_base {
    typedef Eigen::MatrixXd dense_m;
    typedef Eigen::SparseMatrix<double> sparse_m;

    /// get hopping matrix dimension
    int get_msize() const {return m_size_;}; 
    /// get the hopping matrix
    const sparse_m& hopping_m() const { return hopping_m_; }
    /// construct from hopping matrix
    lattice_base(sparse_m in); 
    /// copy constructor 
    lattice_base(lattice_base const& rhs) : hopping_m_(rhs.hopping_m_), m_size_(rhs.m_size_){};
    /// disable moving
    lattice_base(lattice_base && rhs) = delete;

protected:
    /// Hopping matrix
    sparse_m hopping_m_;
    /// Size of the hopping matrix
    size_t m_size_;
};

inline lattice_base::lattice_base(sparse_m in):
   hopping_m_(in),
   m_size_(hopping_m_.rows()) 
{
    if (hopping_m_.rows() != hopping_m_.cols() || hopping_m_.rows() == 0) TRIQS_RUNTIME_ERROR << "Failed to initalize lattice. ";
}

}; // end of namespace FK

//#include "lattice/hypercubic.hpp"
//#include "lattice/triangular.hpp"

