#ifndef __FK_MC_HPP
#define __FK_MC_HPP

#include "common.hpp"
#include "lattice_traits.hpp"
#include "configuration.hpp"
#include <boost/mpi/communicator.hpp>

namespace fk {

// obsevables as a function of MC step - for error estimation
struct observables_t {
    std::vector<double> energies;
};

template <class lattice_t>
class fk_mc 
{
    typedef configuration<lattice_t> config_t;
    boost::mpi::communicator world;
public:
    lattice_t lattice;
    observables_t observables;

    triqs::utility::parameter_defaults solve_defaults() const;

    fk_mc(lattice_t l);
    void solve(utility::parameters p);
};

}; // end of namespace FK

#endif // #ifndef __FK_MC_HPP
