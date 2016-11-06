#ifndef __FK_MC_MEASURE_SPECTRUM_HISTORY_HPP_
#define __FK_MC_MEASURE_SPECTRUM_HISTORY_HPP_

#include "../common.hpp"
#include "../configuration.hpp"
#include <boost/mpi/collectives.hpp>

namespace fk {

struct measure_spectrum_history {
    typedef typename configuration_t::real_array_t  real_array_t;

    configuration_t& config;

    int _Z = 0.0;
    std::vector<std::vector<double>>& _spectrum_history;

    measure_spectrum_history(configuration_t& in, std::vector<std::vector<double>>& spectrum_history);
 
    void measure(double sign);
    void collect_results(boost::mpi::communicator const &c);

};

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURE_SPECTRUM_HISTORY_HPP_
