#include "fk_mc/measures/spectrum.hpp"
#include <boost/mpi/collectives.hpp>

namespace fk {

measure_spectrum::measure_spectrum(configuration_t& in, std::vector<double>& average_spectrum):
config(in), _average_spectrum(average_spectrum) 
{
    _average_spectrum.resize(config.lattice_.get_msize(), 0.0); 
};


void measure_spectrum::accumulate(double sign) 
{
    config.calc_ed(false);
    auto spectrum = config.ed_data_.cached_spectrum;
    for (int i=0; i<_average_spectrum.size(); i++) {
        _average_spectrum[i] = (_average_spectrum[i]*_Z + spectrum[i])/(_Z+1);
        }
    _Z++;
}

void measure_spectrum::collect_results(boost::mpi::communicator const &c)
{
    int sum_Z;
    boost::mpi::reduce(c, _Z, sum_Z, std::plus<int>(), 0);
    std::vector<double> arr(_average_spectrum.size());
    boost::mpi::reduce(c, &_average_spectrum[0], _average_spectrum.size(), &arr[0], std::plus<double>(), 0);
    int n = c.size();
    std::transform(arr.begin(), arr.end(), _average_spectrum.data(), [n](double x){return x / n; });
}


} // end of namespace FK
