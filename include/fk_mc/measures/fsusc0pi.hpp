#ifndef __FK_MC_MEASURE_FSUSC0PI_HPP_
#define __FK_MC_MEASURE_FSUSC0PI_HPP_

#include "../common.hpp"
#include "../configuration.hpp"
#include <boost/mpi/collectives.hpp>

namespace fk {

template <typename lattice_t>
struct measure_nf0pi {
    typedef typename configuration_t::real_array_t  real_array_t;

    const lattice_t& lattice_;
    const configuration_t& config_;

    int _Z = 0.0;
    double average_nf0_ = 0.0;
    double average_nfpi_ = 0.0;
    std::vector<double>& n0_; // n(q=0)
    std::vector<double>& npi_; // n(q=pi)

    measure_nf0pi(const configuration_t& in, const lattice_t& lattice, 
                  std::vector<double>& n0, 
                  std::vector<double>& npi):
        config_(in), lattice_(lattice),
        n0_(n0), npi_(npi) {}
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

};

template <typename lattice_t>
void measure_nf0pi<lattice_t>::accumulate(double sign)
{
    double n0 = config_.f_config_.sum();
    double npi = std::abs(lattice_.FFT_pi(config_.f_config_)); 
    average_nf0_ += n0;
    average_nfpi_ += npi;
    n0_.push_back(n0);
    npi_.push_back(npi);
    _Z++;
}
template <typename lattice_t>
void measure_nf0pi<lattice_t>::collect_results(boost::mpi::communicator const &c)
{
    int sum_Z;
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
        std::cout << "Average nf(q=0): " << nf0_aver / sum_Z << std::endl;
        std::cout << "Average nf(q=pi): " << nfpi_aver / sum_Z << std::endl;
    }
}



} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURE_FSUSC0PI_HPP_
