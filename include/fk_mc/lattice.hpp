#pragma once 

#include "common.hpp"
#include <Eigen/SparseCore>
#include <fftw3.h>

namespace fk {

template <typename T> T myconj(T in);
template <> std::complex<double> myconj<std::complex<melem_type>> (std::complex<melem_type> in) { return std::conj(in); }
template <> double myconj<double> (double in) { return in; }

/** A class to represent a lattice of a finite volume. 
 *  It defines and provides the tight-binding "hopping" matrix in real space.
 */
struct lattice_base {
    typedef Eigen::Matrix<melem_type, Eigen::Dynamic, Eigen::Dynamic> dense_m;
    typedef Eigen::SparseMatrix<melem_type> sparse_m;

    /// get hopping matrix dimension
    size_t get_msize() const {return hopping_m_.rows();}; 
    /// get the total number of orbitals
    size_t get_norbs() const { return norbs_; }
    /// get the volume of the lattice
    size_t get_volume() const { return get_msize() / get_norbs(); } 
    /// get the hopping matrix
    const sparse_m& hopping_m() const { return hopping_m_; }
    /// construct from hopping matrix and the number of orbitals
    lattice_base(sparse_m in, size_t norbitals); 
    /// copy constructor 
    lattice_base(lattice_base const& rhs) : hopping_m_(rhs.hopping_m_){}
    /// disable moving
    lattice_base(lattice_base && rhs) = delete;
    /// The nearest neighbor indices to the given index.
    virtual std::vector<size_t> neighbor_index(size_t index) const = 0;

    virtual size_t ndim() const = 0;

    lattice_base& add_hopping(size_t l, size_t r, melem_type v, bool symmetrize = false) { 
        hopping_m_.insert(l,r) = v;
        if (symmetrize) { 
            hopping_m_.insert(r,l) = myconj(v); 
            }
        return (*this);
        }

protected:
    /// Hopping matrix
    sparse_m hopping_m_;
    /// Number of orbitals (hopping matrix size = volume * n_orbitals)
    size_t norbs_; 
};

inline lattice_base::lattice_base(sparse_m in, size_t norbitals = 1):
   hopping_m_(in),
   norbs_(norbitals) 
{
    if (hopping_m_.rows() != hopping_m_.cols() || hopping_m_.rows() == 0) FKMC_ERROR << "Failed to initalize lattice. ";
}

}; // end of namespace FK


