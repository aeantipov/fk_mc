#include <boost/mpi/collectives.hpp>

#include "spectrum_history.hpp"

namespace fk {

measure_spectrum_history::measure_spectrum_history(const configuration_t& in, std::vector<std::vector<double>>& spectrum_history):
    config(in), _spectrum_history(spectrum_history)
{ 
    _spectrum_history.resize(config.lattice.get_msize());
};

void measure_spectrum_history::accumulate (double sign) 
{
    auto spectrum = config.cached_spectrum;
    for (size_t i=0; i<spectrum.size(); ++i) { _spectrum_history[i].push_back(spectrum(i)); };
    _Z++;
}

void measure_spectrum_history::collect_results(boost::mpi::communicator const &c)
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


} // end of namespace FK
