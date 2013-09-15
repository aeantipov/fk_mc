#ifndef __FK_MC_MEASURE_SPECTRUM_HPP_
#define __FK_MC_MEASURE_SPECTRUM_HPP_

#include "common.hpp"
#include <boost/mpi/collectives.hpp>

namespace fk {

template <class config_t>
struct measure_spectrum {
    typedef typename config_t::real_array_t  real_array_t;

    const config_t& config;

    int _Z = 0.0;
    real_array_t _average_spectrum;

    measure_spectrum(const config_t& in, real_array_t average_spectrum):
        config(in), _average_spectrum(average_spectrum)
        { _average_spectrum.resize(config.get_m_size()); _average_spectrum.setZero();};
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

};

template <class config_t>
void measure_spectrum<config_t>::accumulate (double sign) 
{
    auto spectrum = config.cached_spectrum;
    _average_spectrum = (_average_spectrum*_Z + spectrum)/(_Z+1);
    _Z++;
}

template <class config_t>
void measure_spectrum<config_t>::collect_results(boost::mpi::communicator const &c)
{
    int sum_Z;
    boost::mpi::reduce(c, _Z, sum_Z, std::plus<int>(), 0);
    triqs::arrays::array<double, 1> arr(_average_spectrum.size());
    std::copy(_average_spectrum.data(),_average_spectrum.data()+_average_spectrum.size(),arr.begin());
    boost::mpi::reduce(c, arr, arr, std::plus<triqs::arrays::array<double, 1>>(), 0);
    arr/=c.size();
    std::copy(arr.begin(), arr.end(), _average_spectrum.data());
}

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURE_SPECTRUM_HPP_
