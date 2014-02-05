#ifndef __FK_MC_MEASURE_FOCC_HISTORY_HPP_
#define __FK_MC_MEASURE_FOCC_HISTORY_HPP_

#include <boost/mpi/collectives.hpp>

#include "../configuration.hpp"

namespace fk {

struct measure_focc {
    typedef typename configuration_t::real_array_t real_array_t;
    typedef typename configuration_t::int_array_t int_array_t;

    const configuration_t& config;
    std::vector<std::vector<double>>& focc_;

    measure_focc(const configuration_t& in, std::vector<std::vector<double>>& focc): 
        config(in), focc_(focc)
        { focc_.resize(config.lattice.get_msize()); };
    int _Z = 0;
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

};

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURE_FOCC_HISTORY_HPP_
