#ifndef __FK_MC_MEASURE_FSUSC_HPP_
#define __FK_MC_MEASURE_FSUSC_HPP_

#include "../common.hpp"
#include "../lattice.hpp"
#include "../configuration.hpp"
#include <boost/mpi/collectives.hpp>

namespace fk {

template <typename lattice_t>
struct measure_fsusc {
    typedef typename configuration_t::real_array_t  real_array_t;
    typedef Eigen::ArrayXcd vector_type;
    const lattice_t& lattice_;
    constexpr static size_t D_ = lattice_t::Ndim;
    typedef typename lattice_t::BZPoint BZPoint;
    

    double beta;
    const configuration_t& config_;

    std::vector<BZPoint> qpts_;
    std::vector<std::vector<std::complex<double>>>& nq_history_;
    std::vector<std::vector<double>>& fsuscq_history_;

    const size_t nqpts;

    int _Z = 0.0;

    measure_fsusc(
        double beta, const lattice_t& lattice, const configuration_t& config, 
        std::vector<BZPoint> qpts,
        std::vector<std::vector<std::complex<double>>>& nq_history,
        std::vector<std::vector<double>>& fsuscq_history):
        beta(beta), lattice_(lattice), config_(config), qpts_(qpts), nq_history_(nq_history), fsuscq_history_(fsuscq_history),
        nqpts(qpts.size())
        {
            nq_history_.resize(nqpts);
            fsuscq_history.resize(nqpts);
        };
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);
};

template <typename lattice_t> 
void measure_fsusc<lattice_t>::accumulate(double sign)
{
    
    Eigen::ArrayXcd nf_in = config_.f_config.cast<std::complex<double>>();
    Eigen::ArrayXcd nq = lattice_.FFT(nf_in, FFTW_FORWARD);
    for (size_t i=0; i<nqpts; ++i) { 
        BZPoint q = qpts_[i];
        std::complex<double> nq_val = nq(q.ind_);
        MY_DEBUG(size_t(q) << q <<":" << nq_val << " " << std::abs(nq_val * nq_val));
        nq_history_[i].push_back(nq_val); 
        fsuscq_history_[i].push_back(std::abs(nq_val*nq_val)); 
    };
    MY_DEBUG(lattice_.FFT_pi(config_.f_config));
    _Z++;
}

template <typename lattice_t>
void measure_fsusc<lattice_t>::collect_results(boost::mpi::communicator const &c)
{
    int sum_Z;
    boost::mpi::reduce(c, _Z, sum_Z, std::plus<int>(), 0);
    INFO("F-electron susc and occ number : " << sum_Z << " measures.");
    std::vector<std::complex<double>> nq_history;
    std::vector<double> fsuscq_history;
    for (size_t i=0; i<fsuscq_history_.size(); ++i) {
        fsuscq_history.resize(fsuscq_history_[i].size()*c.size());
        nq_history.resize(nq_history_[i].size()*c.size());
        boost::mpi::gather(c, fsuscq_history_[i].data(), fsuscq_history_[i].size(), fsuscq_history, 0);
        boost::mpi::gather(c, nq_history_[i].data(), nq_history_[i].size(), nq_history, 0);
        fsuscq_history_[i].swap(fsuscq_history);
        nq_history_[i].swap(nq_history);
        };
}

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURE_FSUSC_HPP_
