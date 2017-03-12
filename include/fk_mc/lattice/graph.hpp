#pragma once 

#include <Eigen/SparseCore>

#include "common.hpp"
#include "lattice.hpp"
#include <alps/hdf5.hpp>

namespace fk {

class graph_lattice : public abstract_lattice {
protected:
    using abstract_lattice::hopping_m_;
    using abstract_lattice::norbs_;
public:
    typedef std::vector<int> pos_t;
    graph_lattice(std::vector<pos_t> sites, int norbs = 1);
    graph_lattice(sparse_m hopping, std::vector<pos_t> sites, int norbs = 1);
    virtual size_t ndim() const { return (sites_.size()==0)?0:sites_[0].size();}
    virtual std::vector<size_t> neighbor_index(size_t index) const {}
protected:
    void define_mappings_();
    std::vector<pos_t> sites_;
    std::map<pos_t, int> pos_index_map_;
    pos_t dims_;
};

graph_lattice read_lattice_hdf5(alps::hdf5::archive&, std::string section);

}; // end of namespace FK

