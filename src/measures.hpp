#ifndef __FK_MC_MEASURES_HPP_
#define __FK_MC_MEASURES_HPP_

#include "common.hpp"
#include <boost/mpi/collectives.hpp>

namespace fk {

/*
struct all_measures {
    std::vector<double> weights;
    std::vector<double> energies;
}*/

// a measurement: the magnetization
template <class config_t>
struct measure_energy {
    double beta;
    const config_t& config;

    double Z = 0.0;
    double energy = 0.0;

    double offset_energy = 0.0;
    int offset_nf = 0.0;

    measure_energy(double beta,const config_t& in, double offset_energy):beta(beta),config(in),offset_energy(offset_energy){};
 
    void accumulate(double sign);

    void collect_results(boost::mpi::communicator const &c);

};

template <class config_t>
void measure_energy<config_t>::accumulate (double sign) 
{
    auto evals = config.cached_spectrum;
    auto dm = density_matrix_c(beta,evals,offset_energy);
    MY_DEBUG(dm);
    double weight = __prod(dm);
    Z+=weight;

    auto evs = evals / dm;
    double e_val = __sum(evals);//-config.mu_f*config.get_nf();
    energy += e_val*weight;
    MY_DEBUG(Z);
    MY_DEBUG(energy);
}

template <class config_t>
void measure_energy<config_t>::collect_results(boost::mpi::communicator const &c)
{
    double sum_Z, sum_E;
    boost::mpi::reduce(c, Z, sum_Z, std::plus<double>(), 0);
    boost::mpi::reduce(c, energy, sum_E, std::plus<double>(), 0);

    MY_DEBUG(sum_Z<< " " << sum_E);
    if (c.rank() == 0) {
    INFO("Total energy: " << sum_E / sum_Z);
    }
}

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURES_HPP_
