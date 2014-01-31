#ifndef __FK_MC_MEASURE_SPECTRUM_HPP_
#define __FK_MC_MEASURE_SPECTRUM_HPP_

#include <boost/mpi/communicator.hpp>

#include "../common.hpp"
#include "../configuration.hpp"

namespace fk {

struct measure_spectrum {
    typedef typename configuration_t::real_array_t  real_array_t;

    const configuration_t& config;

    int _Z = 0.0;
    real_array_t& _average_spectrum;

    measure_spectrum(const configuration_t& in, real_array_t& average_spectrum);
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

};

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURE_SPECTRUM_HPP_
