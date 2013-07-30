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

// Free functions
inline real_array_t density_matrix_c(double beta, real_array_t spectrum, double offset_energy)
{
    double exp_offset = exp(beta*offset_energy);
    auto F = triqs::arrays::map(std::function<double(double)>( [beta,offset_energy,exp_offset](double E){return exp_offset+exp(-beta*(E-offset_energy))/exp_offset;} ));
    return F(spectrum);
}


// a measurement: the magnetization
template <class config_t>
struct measure_energy {
    double beta;
    const config_t& config;

    int Z = 0.0;
    double energy = 0.0;

    measure_energy(double beta,const config_t& in):beta(beta),config(in){};
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

};

template <class config_t>
void measure_energy<config_t>::accumulate (double sign) 
{
    auto evals = config.cached_spectrum;
    Z++;

    real_array_t e_nf(evals.shape()[0]);
    triqs::clef::placeholder<0> i_;
    e_nf(i_) << evals(i_) / (1.0+exp(beta*(evals(i_))));

    double e_val = __sum(e_nf) - double(config.mu_f)*config.get_nf();
    energy += e_val;
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
    boost::mpi::reduce(c, Z, sum_Z, std::plus<int>(), 0);
    boost::mpi::reduce(c, energy, sum_E, std::plus<double>(), 0);

    MY_DEBUG(sum_Z<< " " << sum_E);
    if (c.rank() == 0) {
    INFO("Total energy: " << sum_E / sum_Z);
    }
}

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURES_HPP_
