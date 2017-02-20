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
    /// get the total number of orbitals
    size_t get_norbs() const { return norbs_; }
    /// get the hopping matrix
    const sparse_m& hopping_m() const { return hopping_m_; }
    /// construct from hopping matrix and the number of orbitals
    lattice_base(sparse_m in, size_t norbitals); 
    /// copy constructor 
    lattice_base(lattice_base const& rhs) : hopping_m_(rhs.hopping_m_), m_size_(rhs.m_size_){};
    /// disable moving
    lattice_base(lattice_base && rhs) = delete;

    /// The nearest neighbor indices to the given index.
    virtual std::vector<size_t> neighbor_index(size_t index) const = 0;

    virtual size_t ndim() const = 0;

protected:
    /// Hopping matrix
    sparse_m hopping_m_;
    /// Size of the hopping matrix
    size_t m_size_;
    /// Number of orbitals (hopping matrix size = volume * n_orbitals)
    size_t norbs_; 
};

inline lattice_base::lattice_base(sparse_m in, size_t norbitals = 1):
   hopping_m_(in),
   m_size_(hopping_m_.rows()),
   norbs_(norbitals) 
{
    if (hopping_m_.rows() != hopping_m_.cols() || hopping_m_.rows() == 0) FKMC_ERROR << "Failed to initalize lattice. ";
}

}; // end of namespace FK


