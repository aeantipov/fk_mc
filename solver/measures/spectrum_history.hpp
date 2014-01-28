#ifndef __FK_MC_MEASURE_SPECTRUM_HISTORY_HPP_
#define __FK_MC_MEASURE_SPECTRUM_HISTORY_HPP_

#include "../common.hpp"
#include <boost/mpi/collectives.hpp>

namespace fk {

template <class config_t>
struct measure_spectrum_history {
    typedef typename config_t::real_array_t  real_array_t;

    const config_t& config;

    int _Z = 0.0;
    std::vector<std::vector<double>>& _spectrum_history;

    measure_spectrum_history(const config_t& in, std::vector<std::vector<double>>& spectrum_history):
        config(in), _spectrum_history(spectrum_history)
        { _spectrum_history.resize(config.get_m_size());};
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

};

template <class config_t>
void measure_spectrum_history<config_t>::accumulate (double sign) 
{
    auto spectrum = config.cached_spectrum;
    for (size_t i=0; i<spectrum.size(); ++i) { _spectrum_history[i].push_back(spectrum(i)); };
    _Z++;
}

template <class config_t>
void measure_spectrum_history<config_t>::collect_results(boost::mpi::communicator const &c)
{
    int sum_Z;
    boost::mpi::reduce(c, _Z, sum_Z, std::plus<int>(), 0);
    std::vector<double> current_history;
    for (size_t i=0; i<_spectrum_history.size(); ++i) {
        current_history.resize(_spectrum_history[i].size()*c.size());
        boost::mpi::gather(c, _spectrum_history[i].data(), _spectrum_history[i].size(), current_history, 0);
        _spectrum_history[i].swap(current_history);
        };
}

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURE_SPECTRUM_HISTORY_HPP_
