#ifndef __FK_MC_HPP
#define __FK_MC_HPP

#include "common.hpp"
#include "lattice.hpp"
#include "configuration.hpp"

#include <boost/mpi/communicator.hpp>
#include <triqs/mc_tools/mc_generic.hpp>

namespace fk {

struct observables_t {
    typedef configuration_t::dense_m dense_m;
    std::vector<double> energies;
    std::vector<double> d2energies;
    std::vector<double> nf0;
    std::vector<double> nfpi;
    std::vector<double> spectrum;
    std::vector<std::vector<double>> spectrum_history;  // L^D x n_measures size
    std::vector<std::vector<double>> ipr_history;       // L^D x n_measures size
    std::vector<double> stiffness; // n_measures size
    std::vector<std::vector<double>> cond_history;
    std::vector<std::vector<double>> focc_history;      // L^D x n_measures size
    std::vector<std::vector<std::complex<double>>> nq_history;       // nqpts x n_measures size
    std::vector<std::vector<double>> fsuscq_history;    // nqpts x n_measures size
    std::vector<dense_m> eigenfunctions_history;

    void merge(observables_t& rhs);

    void reserve(int n); 
    observables_t() = default;
};

template <typename LatticeType>
class fk_mc : public triqs::mc_tools::mc_generic<double> 
{
    typedef triqs::mc_tools::mc_generic<double> base;
    static_assert(!std::is_same<LatticeType,lattice_base>::value,"Can't construct mc for an unspecified lattice");
    boost::mpi::communicator comm;
public:
    typedef configuration_t config_t;
    typedef LatticeType lattice_type;
    parameters_t p;
    lattice_type lattice;
    configuration_t config;
    //mc_tools::mc_generic<double> mc;
    //const lattice_type& lattice;
    observables_t observables;
    //template <typename MeasureType> void add_measure(MeasureType&& in, std::string name);

    static triqs::utility::parameter_defaults solve_defaults();

    fk_mc(lattice_type l, parameters_t p, bool randomize_config = true);
    void solve(std::vector<double> wgrid_conductivity = {0.0});
    parameters_t& parameters() { return p; } 
};


}; // end of namespace FK

#endif // #ifndef __FK_MC_HPP


