#include "fk_mc/lattice/graph.hpp"

namespace fk { 

lattice_graph::lattice_graph(sparse_m hopping, std::vector<pos_t> sites, int norbs):
    abstract_lattice(hopping, norbs),
    sites_(sites)
{
// define all sorts of mappings
}

} // end of namespace fk
