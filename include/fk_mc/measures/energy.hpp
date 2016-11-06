#ifndef __FK_MC_MEASURE_ENERGY_HPP_
#define __FK_MC_MEASURE_ENERGY_HPP_

#include "../common.hpp"
#include "../configuration.hpp"
#include <boost/mpi/collectives.hpp>

namespace fk {

struct measure_energy {
    typedef typename configuration_t::real_array_t  real_array_t;

    double beta;
    configuration_t& config;

    int _Z = 0.0;
    double _average_energy = 0.0;
    double _average_d2energy = 0.0;
    std::vector<double>& _energies;
    std::vector<double>& _d2energies;
    std::vector<double>& _c_energies;

    measure_energy(double beta, configuration_t& in, std::vector<double>& energies, std::vector<double>& d2energies, std::vector<double>& c_energies):
        beta(beta),config(in), _energies(energies),_d2energies(d2energies),_c_energies(c_energies){};
 
    void measure(double sign);
    void collect_results(boost::mpi::communicator const &c);

};

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURE_ENERGY_HPP_
