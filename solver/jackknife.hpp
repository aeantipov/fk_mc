#ifndef __FK_MC_JACKKNIFE_HPP_
#define __FK_MC_JACKKNIFE_HPP_

#include "binning.hpp"

namespace fk {
namespace jackknife { 

using namespace binning;
typedef bin_data_t bin_data_t;


template <typename Functor, typename iter_t, size_t L> 
struct jackknife_adapter {
typedef typename iter_t::value_type value_type;

    static void jack(Functor F, const std::array<std::pair<iter_t,iter_t>, L> &data )
    {
        size_t size = std::distance(data[0].first,data[0].second);
        std::array<value_type, L> x_mean;
        for (size_t i=0; i<L; ++i) x_mean[i] = std::get<_MEAN>(calc_stats(data[i].first, data[i].second));
        double U_0 = triqs::tuple::apply(F,x_mean); 
        MY_DEBUG("U0: " << U_0);

        std::vector<double> U(size);
        std::array<iter_t, L> iters; 
        for (size_t i=0; i<L; ++i) { iters[i]=data[i].first; };

        std::array<double, L> xy_means;
        for (size_t j=0; j<size; ++j) { 
        xy_means.fill(0.0);
        for (size_t i=0; i<L; ++i) { 
            xy_means[i] = std::accumulate(data[i].first, iters[i], 0.0, std::plus<double>());
            iters[i]++;
            xy_means[i]+= std::accumulate(iters[i], data[i].second, 0.0, std::plus<double>());
            xy_means[i]/=(size-1);
        }
        double U_j = triqs::tuple::apply(F,xy_means);
        U[j] = U_j;
        //MY_DEBUG("U"<<j+1<<": " << U_j);
        }
        auto U_stats = calc_stats(U.begin(),U.end());
        double U_bar = std::get<_MEAN>(U_stats);
        double U_average = U_0 - (size-1)*(U_bar - U_0);
        double dU = (size-1)*std::get<_SQERROR>(U_stats);
        MY_DEBUG(U_0 << " " << U_average << "+/-" << dU);
        //triqs::array<value_type, 1> jack_means;
    }

};

} // end of namespace jackknife
} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_JACKKNIFE_HPP_

