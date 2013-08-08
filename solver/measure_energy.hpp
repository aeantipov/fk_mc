#ifndef __FK_MC_MEASURES_HPP_
#define __FK_MC_MEASURES_HPP_

#include "common.hpp"
#include <boost/mpi/collectives.hpp>

namespace fk {

template <class config_t>
struct measure_energy {
    double beta;
    const config_t& config;

    int _Z = 0.0;
    double _average_energy = 0.0;
    std::vector<double>& _energies;

    measure_energy(double beta,const config_t& in, std::vector<double>& energies):beta(beta),config(in), _energies(energies){};
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

};

template <class config_t>
void measure_energy<config_t>::accumulate (double sign) 
{
    auto spectrum = config.cached_spectrum;
    _Z++;

    real_array_t e_nf(spectrum.shape()[0]);
    triqs::clef::placeholder<0> i_;
    e_nf(i_) << spectrum(i_) / (1.0+exp(beta*(spectrum(i_))));

    double e_val = sum(e_nf) - double(config.mu_f)*config.get_nf();
    _average_energy += e_val;
    _energies.push_back(e_val);
    /*MY_DEBUG(__sum(e_nf));
    MY_DEBUG(e_val);
    MY_DEBUG(Z);
    MY_DEBUG(energy);
    */
}

template <class config_t>
void measure_energy<config_t>::collect_results(boost::mpi::communicator const &c)
{
    int sum_Z;
    double sum_E;
    boost::mpi::reduce(c, _Z, sum_Z, std::plus<int>(), 0);
    boost::mpi::reduce(c, _average_energy, sum_E, std::plus<double>(), 0);

    std::vector<double> energies(_energies.size()*c.size());
    boost::mpi::gather(c, _energies.data(), _energies.size(), energies, 0);
    _energies.swap(energies);

    if (c.rank() == 0) {
    INFO("Total energy: " << sum_E / sum_Z);
    }
}

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURES_HPP_
