#include "fk_mc/lattice/graph.hpp"

namespace fk { 

template <typename T>
std::ostream& pprint(std::ostream& out, T x, std::string split_s = " ") { for (auto x1 : x) out << x1 << split_s; out << std::endl;}

graph_lattice::graph_lattice(std::vector<pos_t> sites, int norbs):
    abstract_lattice(norbs),
    sites_(sites)
{
    define_mappings_();
}

graph_lattice::graph_lattice(sparse_m hopping, std::vector<pos_t> sites, int norbs):
    abstract_lattice(hopping, norbs),
    sites_(sites)
{
    define_mappings_();
}

void graph_lattice::define_mappings_()
{
    int p = 0;
    for (pos_t pos : sites_) {
        pos_index_map_[pos] = p;
        ++p;
        assert(pos.size() == ndim());
    }
}


graph_lattice read_lattice_hdf5(alps::hdf5::archive& ar, std::string section)
{
    std::cout << "Loading /" << section << std::endl;
    std::vector<int> dims;
    ar[section + "/dims"] >> dims;
    int norbs; 
    ar[section + "/norbitals"] >> norbs;
    int ndims = dims.size();
    int npts = std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<int>());
    std::cout << ndims << " dimensions, [" << norbs << " orbitals, " << npts << " points." << std::endl;
    int m_size = npts * norbs;

    std::vector<graph_lattice::pos_t> sites;
    ar[section + "/coordinates"] >> sites;
    for (auto x : sites) pprint(std::cout, x); 

    std::vector<std::vector<melem_type>> m_in_data;
    ar[section + "/hopping"] >> m_in_data;
    typename graph_lattice::sparse_m hopping_m(m_size, m_size);
    for (auto v : m_in_data) { 
        int p1 = std::real(v[0]); 
        int p2 = std::real(v[1]);
        melem_type m = v[2]; 

        hopping_m.coeffRef(p1, p2) = m;
    }   

    graph_lattice l(hopping_m, sites, norbs);
    return l;
}

} // end of namespace fk
