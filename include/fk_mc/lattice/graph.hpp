#pragma once 

#include <Eigen/SparseCore>

#include "common.hpp"
#include "lattice.hpp"

namespace fk {

class lattice_graph : public abstract_lattice {
protected:
    using abstract_lattice::hopping_m_;
    using abstract_lattice::norbs_;
public:
    typedef std::vector<size_t> pos_t;
    virtual size_t ndims() { return ndims_ ;}
    lattice_graph(sparse_m hopping, std::vector<pos_t> sites, int norbs = 1);// : base(hopping, norbs), sites_(sit)
protected:
    int ndims_;
    std::vector<pos_t> sites_;
    pos_t dims_;
};

}; // end of namespace FK

