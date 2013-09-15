#ifndef __FK_MC_MEASURE_D2E_HPP_
#define __FK_MC_MEASURE_D2E_HPP_

#include "common.hpp"
#include <boost/mpi/collectives.hpp>

namespace fk {

template <class config_t>
struct measure_d2energy {
    typedef typename config_t::real_array_t  real_array_t;

    double beta;
    const config_t& config;

    int _Z = 0.0;
    double _average_d2energy = 0.0;
    std::vector<double>& _d2energies;

    measure_d2energy(double beta,const config_t& in, std::vector<double>& d2energies):beta(beta),config(in), _d2energies(d2energies){};
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

};

template <class config_t>
void measure_d2energy<config_t>::accumulate (double sign) 
{
    auto spectrum = config.cached_spectrum;
    _Z++;

    double d2e_val = d2e_nf.sum()/2.0;
    _average_d2energy += d2e_val;
    _d2energies.push_back(d2e_val);
}

template <class config_t>
void measure_d2energy<config_t>::collect_results(boost::mpi::communicator const &c)
{
    int sum_Z;
    double sum_d2E;
    boost::mpi::reduce(c, _Z, sum_Z, std::plus<int>(), 0);
    boost::mpi::reduce(c, _average_d2energy, sum_d2E, std::plus<double>(), 0);

    std::vector<double> d2energies(_d2energies.size()*c.size());
    boost::mpi::gather(c, _d2energies.data(), _d2energies.size(), d2energies, 0);
    _d2energies.swap(d2energies);

    if (c.rank() == 0) {
    INFO("Total d2energy: " << sum_d2E / sum_Z);
    }
}

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURE_D2E_HPP_
