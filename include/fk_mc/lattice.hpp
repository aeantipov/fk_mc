#pragma once 

#include "common.hpp"
#include <Eigen/SparseCore>
#include <fftw3.h>

namespace fk {

/** A class to represent a lattice of a finite volume. 
 *  It defines and provides the tight-binding "hopping" matrix in real space.
 */
struct abstract_lattice {
    typedef Eigen::Matrix<melem_type, Eigen::Dynamic, Eigen::Dynamic> dense_m;
    typedef Eigen::SparseMatrix<melem_type> sparse_m;

    /// get hopping matrix dimension
    size_t msize() const {return hopping_m_.rows();};
    /// get the total number of orbitals
    size_t norbs() const { return norbs_; }
    /// get the volume of the lattice
    size_t volume() const { return size_t(msize() / norbs()); }
    /// get the hopping matrix
    const sparse_m& hopping_m() const { return hopping_m_; }
    /// construct from hopping matrix and the number of orbitals
    abstract_lattice(sparse_m in, size_t norbitals);
    /// copy constructor 
    abstract_lattice(abstract_lattice const& rhs) : hopping_m_(rhs.hopping_m_), norbs_(rhs.norbs_){}
    /// disable moving
    abstract_lattice(abstract_lattice && rhs) = delete;
    /// The nearest neighbor indices to the given index.
    virtual std::vector<size_t> neighbor_index(size_t index) const = 0;

    virtual size_t ndim() const = 0;

    abstract_lattice& add_hopping(size_t l, size_t r, melem_type v, bool symmetrize = false);

protected:
    /// Hopping matrix
    sparse_m hopping_m_;
    /// Number of orbitals (hopping matrix size = volume * n_orbitals)
    size_t norbs_; 
};

inline abstract_lattice::abstract_lattice(sparse_m in, size_t norbitals = 1):
   hopping_m_(in),
   norbs_(norbitals) 
{
    if (hopping_m_.rows() != hopping_m_.cols() || hopping_m_.rows() == 0) FKMC_ERROR << "Failed to initialize lattice. ";
}

}; // end of namespace FK


