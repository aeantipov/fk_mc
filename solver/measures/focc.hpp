#ifndef __FK_MC_MEASURE_FSUSC_HPP_
#define __FK_MC_MEASURE_FSUSC_HPP_

#include "../common.hpp"
#include <boost/mpi/collectives.hpp>

namespace fk {

template <class config_t>
struct measure_focc {
    typedef typename config_t::real_array_t real_array_t;
    typedef typename config_t::int_array_t int_array_t;

    const config_t& config;
    std::vector<std::vector<double>>& focc_;

    measure_focc(const config_t& in, std::vector<std::vector<double>>& focc): 
        config(in), focc_(focc)
        { focc_.resize(config.get_m_size()); };
    int _Z = 0;
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

};

template <class config_t>
void measure_focc<config_t>::accumulate (double sign) 
{
    const auto& focc_current = config.f_config;
    for (size_t i=0; i<focc_current.size(); ++i) { focc_[i].push_back(focc_current(i)); };
    _Z++;
}

template <class config_t>
void measure_focc<config_t>::collect_results(boost::mpi::communicator const &c)
{
    int sum_Z;
    boost::mpi::reduce(c, _Z, sum_Z, std::plus<int>(), 0);
    std::vector<double> focc_tmp;
    for (size_t i=0; i<focc_.size(); ++i) {
        focc_tmp.resize(focc_[i].size()*c.size());
        boost::mpi::gather(c, focc_[i].data(), focc_[i].size(), focc_tmp, 0);
        focc_[i].swap(focc_tmp);
        };
}

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURE_FSUSC_HPP_
