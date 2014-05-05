#ifndef __FK_MC_MEASURE_IPR_HPP__
#define __FK_MC_MEASURE_IPR_HPP__

#include "../common.hpp"
#include "../configuration.hpp"
#include <boost/mpi/collectives.hpp>
#include <Eigen/Dense>

namespace fk {

template <typename lattice_t>
struct measure_ipr {
    typedef typename configuration_t::real_array_t  real_array_t;

    const lattice_t& lattice_;
    configuration_t& config_;

    int _Z = 0.0;
    std::vector<std::vector<double>>& ipr_vals_;

    measure_ipr(configuration_t& in, const lattice_t& lattice, 
                  std::vector<std::vector<double>>& ipr_vals); 
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

};


template <typename lattice_t>
measure_ipr<lattice_t>::measure_ipr(configuration_t& in, const lattice_t& lattice, std::vector<std::vector<double>>& ipr_vals): 
    config_(in), 
    lattice_(lattice),
    ipr_vals_(ipr_vals)
{
    ipr_vals_.resize(config_.lattice_.get_msize());
}

template <typename lattice_t>
void measure_ipr<lattice_t>::accumulate(double sign)
{
    if (int(config_.ed_data_.status) < (ed_cache::full)) throw std::logic_error("Need calculated eigenvalues and eigenvectors");  

    const auto& evecs = config_.ed_data_.cached_evecs;
    const auto& evals = config_.ed_data_.cached_spectrum;

    Eigen::VectorXd ipr(evals.size());
    auto rphi2 = evecs.colwise().squaredNorm();
    for (size_t i=0; i<evals.size(); i++) {
        const auto &t = evecs.col(i);
        auto rphi4_val = t.template lpNorm<4>();
        ipr(i) = rphi4_val/rphi2(i);
    }
    for (size_t i=0; i<evals.size(); ++i) { ipr_vals_[i].push_back(ipr(i)); };
    _Z++;
}
template <typename lattice_t>
void measure_ipr<lattice_t>::collect_results(boost::mpi::communicator const &c)
{
    int sum_Z;
    boost::mpi::reduce(c, _Z, sum_Z, std::plus<int>(), 0);
    std::vector<double> current_ipr;
    for (size_t i=0; i<ipr_vals_.size(); ++i) {
        current_ipr.resize(ipr_vals_[i].size()*c.size());
        boost::mpi::gather(c, ipr_vals_[i].data(), ipr_vals_[i].size(), current_ipr, 0);
        ipr_vals_[i].swap(current_ipr);
        };
}

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURE_IPR_HPP__
