#ifndef __FK_MC_JACKKNIFE_HPP_
#define __FK_MC_JACKKNIFE_HPP_

#include "binning.hpp"

namespace fk {
namespace jackknife { 

#define BINNING_RANGE ((15)(14)(13)(12)(11)(10)(9)(8)(7)(6)(5)(4)(3)(2)(1)(0))

using namespace binning;
typedef bin_data_t bin_data_t;

typedef std::tuple<size_t, double, double, double> jack_stats_t;

struct jackknife_adapter
{
//typedef typename iter_t::value_type value_type;
    template <typename iter_pair>
    static double _calc_mean_from_pair(const iter_pair& in)
    {
        return std::accumulate(in.first, in.second,0.0)/std::distance(in.first,in.second);
    } 

    template <typename Functor, typename iter_t, size_t L>
    static bin_stats_t jack(Functor F, const std::array<std::pair<iter_t,iter_t>, L> &data )
    {
        typedef std::pair<iter_t,iter_t> iter_pair_t;

        size_t size = std::distance(std::get<0>(data).first,std::get<0>(data).second);
        std::array<double, L> x_mean;
        for (size_t i=0; i<L; ++i) x_mean[i] = std::get<_MEAN>(calc_stats(data[i].first, data[i].second));
        double U_0 = triqs::tuple::apply(F,x_mean); 
        //MY_DEBUG("U0: " << U_0);

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
        }
        auto U_stats = calc_stats(U.begin(),U.end());
        double U_bar = std::get<_MEAN>(U_stats);
        double U_average = U_0 - (size-1)*(U_bar - U_0);
        double dU = (size-1)*std::get<_SQERROR>(U_stats);
        double var_U = dU*dU*size; 
        MY_DEBUG(U_0 << " " << U_average << "+/-" << dU);
        return std::make_tuple(size, U_average, var_U, dU);
    }

};

template <class Functor, size_t L, size_t total_bins_left, size_t current_bin = 0>
struct jackknife_accumulator {
    template <typename iter_t>
    static bin_data_t accumulate_jackknife(Functor f, const std::array<std::pair<iter_t,iter_t>, L> &data) {
        typedef binning::binned_iterator<iter_t,current_bin> bin_it;
        std::array<std::pair<bin_it, bin_it>, L> range;
        for (size_t i=0; i<L; i++) range[i]=find_bin_range(data[i].first, data[i].second, bin_it());
        auto stats = jackknife_adapter::jack<Functor,bin_it,L>(f,range);
        auto next  = jackknife_accumulator<Functor, L, total_bins_left - 1, current_bin + 1>::accumulate_jackknife(f,data);
        next.insert(next.begin(),stats);
        return next;
    };
};

template <class Functor, size_t L, size_t current_bin>
struct jackknife_accumulator<Functor,L,0,current_bin> {
    template <typename iter_t>
    static bin_data_t accumulate_jackknife(Functor f, const std::array<std::pair<iter_t,iter_t>, L> &data) {
        bin_data_t out;
        out.reserve(current_bin);
        typedef binning::binned_iterator<iter_t,current_bin> bin_it;
        std::array<std::pair<bin_it, bin_it>, L> range;
        for (size_t i=0; i<L; i++) range[i]=find_bin_range(data[i].first, data[i].second, bin_it());
        auto stats = jackknife_adapter::jack<Functor,bin_it,L>(f,range);
        out.push_back(stats);
        return out;
    };
};

template <typename Functor, size_t L, typename iter_t>
bin_data_t accumulate_jackknife(Functor F, const std::array<std::pair<iter_t,iter_t>, L>& in, size_t bin_depth) {
    #define MACRO(r, p) \
    if (BOOST_PP_SEQ_ELEM(0, p) == bin_depth) \
        return jackknife_accumulator<Functor,L,BOOST_PP_SEQ_ELEM(0, p)>::template accumulate_jackknife<iter_t>(F,in);
    BOOST_PP_SEQ_FOR_EACH_PRODUCT(MACRO, BINNING_RANGE)
    #undef MACRO
    TRIQS_RUNTIME_ERROR << "bin_depth =" << bin_depth << "> compiled bin size";
    return {std::make_tuple(0, std::nan(""),std::nan(""),std::nan(""))};
}

/*
template <typename iter_t, typename Ret, typename ... Args>
bin_data_t accumulate_jackknife(std::function<Ret(Args...)> F, const std::array<std::pair<iter_t,iter_t>, sizeof...(Args)>& in, size_t bin_depth) {
    typedef std::function<Ret(Args...)> Functor;
    constexpr size_t L = sizeof...(Args);
    return accumulate_jackknife<Functor, L, iter_t>(std::forward<Functor>(F),in,bin_depth);
    }
*/

template <typename iter_t, typename Ret, typename ... Args>
bin_data_t accumulate_jackknife(std::function<Ret(Args...)> F, std::initializer_list<std::pair<iter_t,iter_t>> in, size_t bin_depth) {
    typedef std::function<Ret(Args...)> Functor;
    constexpr size_t L = sizeof...(Args);
    std::array<std::pair<iter_t,iter_t>, L> d1(in);

    return accumulate_jackknife<Functor, L, iter_t>(std::forward<Functor>(F),d1,bin_depth);
    }






template <typename Functor, size_t L, typename container_t>
bin_data_t accumulate_jackknife(Functor &&F, const std::array<container_t, L>& in, size_t bin_depth) {
    typedef typename std::add_const<container_t>::type::const_iterator iter_t;
    std::array<std::pair<iter_t,iter_t>, L> ranges;
    for (size_t i=0; i<L; ++i) ranges[i]=std::make_pair(in[i].begin(),in[i].end());
    return accumulate_jackknife(F,ranges,bin_depth);
}



} // end of namespace jackknife
} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_JACKKNIFE_HPP_

