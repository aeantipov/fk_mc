#ifndef __FK_MC_MEASURE_SPECTRUM_HPP_
#define __FK_MC_MEASURE_SPECTRUM_HPP_

#include <boost/mpi/communicator.hpp>

#include "../common.hpp"
#include "../configuration.hpp"

namespace fk {

struct measure_spectrum {
    configuration_t& config;

    int _Z = 0.0;
    std::vector<double>& _average_spectrum;

    measure_spectrum(configuration_t& in, std::vector<double>& average_spectrum);
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

};

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURE_SPECTRUM_HPP_
