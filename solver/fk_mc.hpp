#ifndef __FK_MC_HPP
#define __FK_MC_HPP

#include "common.hpp"
#include "lattice.hpp"
#include "configuration.hpp"

#include <boost/mpi/communicator.hpp>
#include <triqs/mc_tools/mc_generic.hpp>

namespace fk {

struct observables_t {
    std::vector<double> energies;
    std::vector<double> d2energies;
    std::vector<double> nf0;
    std::vector<double> nfpi;
    typename configuration_t::real_array_t spectrum;
    std::vector<std::vector<double>> spectrum_history; // L^D x n_measures size
    std::vector<std::vector<double>> focc_history;     // L^D x n_measures size
    std::vector<std::vector<std::complex<double>>> nq_history;       // nqpts x n_measures size
    std::vector<std::vector<double>> fsuscq_history;   // nqpts x n_measures size
};

class fk_mc 
{
    typedef configuration_t config_t;
    boost::mpi::communicator world;
public:
    utility::parameters p;
    configuration_t config;
    mc_tools::mc_generic<double> mc;
    const lattice_base& lattice;
    observables_t observables;
    template <typename MeasureType> void add_measure(MeasureType&& in, std::string name);

    static triqs::utility::parameter_defaults solve_defaults();

    fk_mc(const lattice_base& l, utility::parameters p);
    void solve();
};

template <typename MeasureType>
void fk_mc::add_measure(MeasureType&& in, std::string name)
{
    mc.add_measure(in,name);
}

}; // end of namespace FK

#endif // #ifndef __FK_MC_HPP
