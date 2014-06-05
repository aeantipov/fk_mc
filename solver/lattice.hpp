#pragma once 

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

inline lattice_base::lattice_base(sparse_m &&in):
   hopping_m(in),
   m_size_(hopping_m.rows()) 
{
    if (hopping_m.rows() != hopping_m.cols() || hopping_m.rows() == 0) TRIQS_RUNTIME_ERROR << "Failed to initalize lattice. ";
}

}; // end of namespace FK

//#include "lattice/hypercubic.hpp"
//#include "lattice/triangular.hpp"

