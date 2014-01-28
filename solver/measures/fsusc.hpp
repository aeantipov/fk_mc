#ifndef __FK_MC_MEASURE_FSUSC_HPP_
#define __FK_MC_MEASURE_FSUSC_HPP_

#include "../common.hpp"
#include <boost/mpi/collectives.hpp>

namespace fk {

template <class config_t>
struct measure_fsusc {
    typedef typename config_t::real_array_t real_array_t;
    typedef typename config_t::int_array_t int_array_t;
    typedef std::array<double,config_t::lattice_t::Ndim> BZPoint;

    template <typename VType> static VType shift(const VType& in, int n);

    std::vector<BZPoint> qpts;
    const config_t& config;

    int _Z = 0;
    std::vector<std::vector<double>>& _average_fsusc;

    measure_fsusc(const config_t& in, std::vector<std::vector<double>>& average_fsusc, std::vector<BZPoint> qpts):
        config(in), _average_fsusc(average_fsusc),qpts(qpts)
        { average_fsusc.resize(qpts.size()); } ;
 
    void accumulate(double sign);
    void collect_results(boost::mpi::communicator const &c);

};

template <class config_t>
template <typename VType>
VType measure_fsusc<config_t>::shift(const VType& in, int n)
{
    int size = in.size();
    n = (n+size)%size;
    VType out(size);
    out.head(n) = in.tail(n);
    out.tail(size-n) = in.head(size-n);
    return out;
}

template <class config_t>
void measure_fsusc<config_t>::accumulate (double sign) 
{
    _Z++;
    real_array_t f_config = config.f_config.template cast<double>();
    real_array_t f_config_shift = f_config;
    real_array_t ff_cor (f_config.size());
    real_array_t ff_cor2 (f_config.size());
    for (size_t i=0; i<f_config.size(); i++) { 
        MY_DEBUG(f_config_shift.transpose());
        auto ff = f_config*f_config_shift;
        ff_cor(i) = ff.mean();
        ff_cor2(i) = f_config(0)*f_config(i);
        f_config_shift = shift(f_config_shift,1);
        };
    MY_DEBUG(ff_cor.transpose());
    MY_DEBUG(ff_cor2.transpose());
}

template <class config_t>
void measure_fsusc<config_t>::collect_results(boost::mpi::communicator const &c)
{
}

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_MEASURE_FSUSC_HPP_
