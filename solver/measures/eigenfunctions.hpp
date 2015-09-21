#pragma once

#include "../common.hpp"
#include "../configuration.hpp"
#include <boost/mpi/collectives.hpp>

namespace fk {

struct measure_eigenfunctions {
    typedef typename configuration_t::real_array_t  real_array_t;
    typedef typename configuration_t::dense_m dense_m;

    configuration_t& config;

    int _Z = 0.0;
    std::vector<dense_m>& eigenfunctions_;

    measure_eigenfunctions(configuration_t& in, std::vector<dense_m>& eigenfunctions);
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

};

} // end of namespace fk

