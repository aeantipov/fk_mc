#ifndef __FK_MC_HPP
#define __FK_MC_HPP

#include "common.hpp"
#include "lattice.hpp"
#include "configuration.hpp"
#include <boost/mpi/communicator.hpp>

namespace fk {

struct observables_t {
    std::vector<double> energies;
    std::vector<double> d2energies;
    typename configuration_t::real_array_t spectrum;
    std::vector<std::vector<double>> spectrum_history; // n_eigenvalues x n_measures size
    std::vector<std::vector<double>> focc_history;     // n_focc x n_measures size
};

class fk_mc 
{
    typedef configuration_t config_t;
    boost::mpi::communicator world;
public:
    const lattice_base& lattice;
    observables_t observables;

    triqs::utility::parameter_defaults solve_defaults() const;

    fk_mc(lattice_base l);
    void solve(utility::parameters p);
};

}; // end of namespace FK

#endif // #ifndef __FK_MC_HPP
