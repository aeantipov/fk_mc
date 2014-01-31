#include "spectrum.hpp"
#include <boost/mpi/collectives.hpp>

namespace fk {

measure_spectrum::measure_spectrum(const configuration_t& in, real_array_t& average_spectrum):
config(in), _average_spectrum(average_spectrum) 
{
    _average_spectrum.resize(config.lattice.get_msize()); 
    _average_spectrum.setZero();
};


void measure_spectrum::accumulate (double sign) 
{
    auto spectrum = config.cached_spectrum;
    _average_spectrum = (_average_spectrum*_Z + spectrum)/(_Z+1);
    _Z++;
}

void measure_spectrum::collect_results(boost::mpi::communicator const &c)
{
    int sum_Z;
    boost::mpi::reduce(c, _Z, sum_Z, std::plus<int>(), 0);
    triqs::arrays::array<double, 1> arr(_average_spectrum.size());
    std::copy(_average_spectrum.data(),_average_spectrum.data()+_average_spectrum.size(),arr.begin());
    boost::mpi::reduce(c, arr, arr, std::plus<triqs::arrays::array<double, 1>>(), 0);
    arr/=c.size();
    std::copy(arr.begin(), arr.end(), _average_spectrum.data());
}


} // end of namespace FK
