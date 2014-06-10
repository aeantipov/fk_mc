#ifndef __FK_MC_MEASURE_FSUSC0PI_HPP_
#define __FK_MC_MEASURE_FSUSC0PI_HPP_

#include "../common.hpp"
#include "../configuration.hpp"
#include <boost/mpi/collectives.hpp>

namespace fk {

template <typename Lattice>
struct measure_polarization {
    typedef typename configuration_t::real_array_t  real_array_t;

    const Lattice& lattice_;
    configuration_t& config_;

    std::vector<std::complex<double>> pol_v;
    Eigen::VectorXcd phase_v;

    measure_polarization(configuration_t& in, const Lattice& lattice, int dir = 0): 
        config_(in), lattice_(lattice){
            phase_v.resize(lattice.get_msize());
            for (int j=0; j<lattice.get_msize(); j++) { 
                phase_v[j] = exp(I*double(lattice_.index_to_pos(j)[dir])*2.*M_PI / double(lattice_.dims[dir]));
                }
        }
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

};

template <typename Lattice>
void measure_polarization<Lattice>::accumulate(double sign)
{
    config_.calc_ed(true);

    const auto& evecs = config_.ed_data().cached_evecs;
    const auto& evals = config_.ed_data().cached_spectrum;

    std::complex<double> pol = 0.0;
    double nc = 0.0;

    for (size_t i=0; i<evals.size() && config_.ed_data().cached_weights[i]>=1e-5; i++) {
        const auto &ev = evecs.col(i);
        std::complex<double> x = (ev.array()*ev.array()).matrix().transpose()*phase_v;
        FKDEBUG(i << " " << config_.ed_data().cached_weights[i] << " : " << x << " -> " << x*config_.ed_data().cached_weights[i]);
        FKDEBUG(ev.transpose() << std::endl);
        //FKDEBUG(phase_v.transpose());
        pol+=x*config_.ed_data().cached_weights[i];
    }
    FKDEBUG(pol);
    exit(0);

    pol_v.push_back(pol);
    //_Z++;
}

template <typename Lattice>
void measure_polarization<Lattice>::collect_results(boost::mpi::communicator const &c)
{
   // int sum_Z;
/*
    double nf0_aver = 0, nfpi_aver = 0;
    boost::mpi::reduce(c, _Z, sum_Z, std::plus<int>(), 0);
    boost::mpi::reduce(c, average_nf0_, nf0_aver, std::plus<double>(), 0);
    boost::mpi::reduce(c, average_nfpi_, nfpi_aver, std::plus<double>(), 0);

    size_t datasize = n0_.size()*c.size();
    std::vector<double> tmp(datasize);

    boost::mpi::gather(c, n0_.data(), n0_.size(), tmp, 0);
    n0_.swap(tmp);
    c.barrier();

    tmp.resize(datasize);
    boost::mpi::gather(c, npi_.data(), npi_.size(), tmp, 0);
    npi_.swap(tmp);
    c.barrier();
    
    if (c.rank() == 0) {
    INFO("Average nf(q=0): " << nf0_aver / sum_Z);
    INFO("Average nf(q=pi): " << nfpi_aver / sum_Z);
    }
*/
}

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURE_FSUSC0PI_HPP_
