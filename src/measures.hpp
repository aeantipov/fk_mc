#ifndef __FK_MC_MEASURES_HPP_
#define __FK_MC_MEASURES_HPP_

#include "common.hpp"


namespace fk {

/*
struct all_measures {
    std::vector<double> weights;
    std::vector<double> energies;
}*/

// a measurement: the magnetization
template <class lattice>
struct measure_energy {

    typedef configuration<lattice> config_t;
    double beta;
    const config_t& config;

    double z = 0.0;
    double e = 0.0;

    double offset_energy = 0.0;

    measure_energy(double beta,const config_t& in):beta(beta),config(in){};
 
    void accumulate(double sign);

    void collect_results(boost::mpi::communicator const &c){};

//        double sum_Z, sum_M;
//    boost::mpi::reduce(c, Z, sum_Z, std::plus<double>(), 0);
//    boost::mpi::reduce(c, M, sum_M, std::plus<double>(), 0);

//        if (c.rank() == 0) {
//      std::cout << "Magnetization: " << sum_M / sum_Z << std::endl << std::endl;

};

template <class lattice>
void measure_energy<lattice>::accumulate (double sign) 
{
    DEBUG(sign);
    z+=sign;

    auto evals = config.cached_spectrum;
    double ground_energy = evals(0);
    auto F = triqs::arrays::map(std::function<double(double)>( [this,ground_energy](double E){return ((E-ground_energy)/(1.0+exp(-beta*(E-ground_energy))));}));
    evals = F(evals); 
    double e = __sum(evals)-config.mu_f*config.get_nf();
//    e+=
}

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURES_HPP_
