#ifndef __FK_MC_BINNING_HPP_
#define __FK_MC_BINNING_HPP_

#include "common.hpp"
#include "triqs_extra.hpp"

#include <cmath>
#include <boost/iterator/iterator_facade.hpp>
#include <triqs/utility/exceptions.hpp>
#include <boost/preprocessor.hpp>
#include <iterator>

namespace fk {

/** Binned iterator takes bin as an input parameter. If bin = 0 -> normal iterator, if bin>0 dereference returns averaged over (2^bin) values. */
template <typename iter_t, size_t D>
struct binned_const_iterator : 
    public boost::iterator_facade<binned_const_iterator<iter_t,D>, const double, boost::forward_traversal_tag, double>
{
    constexpr static size_t step = binned_const_iterator<iter_t,D-1>::step*2;
    iter_t p;
    binned_const_iterator(iter_t start):p(start){};
    constexpr double dereference() const { return (*binned_const_iterator<iter_t,D-1>(p) + *binned_const_iterator<iter_t,D-1>(p+binned_const_iterator<iter_t,D-1>::step))/2.; };
    bool equal(binned_const_iterator const &rhs) const { return (p == rhs.p); };
    void increment() { std::advance(p,step); };
};

template <typename iter_t>
struct binned_const_iterator<iter_t,0>: public boost::iterator_facade<binned_const_iterator<iter_t,0>, double, boost::forward_traversal_tag, double>
{
    constexpr static size_t step = 1;
    iter_t p;
    binned_const_iterator(iter_t start):p(start){};
    bool equal(binned_const_iterator const &rhs) const { return (p == rhs.p); };
    constexpr double dereference() const { return *p; };
    void increment() { p++; };
};

template <typename container_t, size_t maxD = 15>
struct binning_adapter {
    //typedef typename container_t::iterator iter_t;
    template <size_t bin, typename iter_t> static std::tuple<double, double, double> binning(iter_t begin, iter_t end, std::random_access_iterator_tag) {
        if (bin > maxD) TRIQS_RUNTIME_ERROR << "Can't bin " << bin << " times, max = " << maxD;
        int size = std::distance(begin, end);
        int step = binned_const_iterator<iter_t,bin>::step;
        if (step > size) TRIQS_RUNTIME_ERROR << "Can't bin with depth = " << bin << ", binning step(" << step << ")> container size (" << size << ")";
        int nsteps = size/step;
        binned_const_iterator<iter_t,bin> binned_begin(begin), binned_end(begin); std::advance(binned_end,nsteps);
        return calc_stats(binned_begin, binned_end, nsteps);
        /*
        double mean, variance, stddev;
        std::tie(mean,variance,stddev) = calc_stats(
        MY_DEBUG(bin << " " << mean << " " << std::sqrt(variance/nsteps));
        return std::make_tuple(mean, variance, std::sqrt(variance/nsteps));
        */
        };

/*
    #define RANGE ((0)(1)(2)(3)(4)(5)(6)(7)(8)(9)(10)(11)(12)(13))
    #define MACRO(r, p) \
    if (BOOST_PP_SEQ_ELEM(0, p) == bin) \
        auto it = binned_const_iterator<BOOST_PP_SEQ_ELEM(0, p)>(in.begin());
    BOOST_PP_SEQ_FOR_EACH_PRODUCT(MACRO, RANGE)
    #undef MACRO
    #undef RANGE

//       binned_const_iterator<1> begin(&*in.begin());
*/
    template <typename iter_t>
    static std::tuple<double, double, double> calc_stats(const iter_t& begin, const iter_t& end, size_t size) {
        double mean = std::accumulate(begin,end,0.0, std::plus<double>())/size; 
        double variance = std::accumulate(begin,end,0.0, [mean](double x, double y){return x + (y-mean)*(y-mean);})/(size-1);
        INFO_NONEWLINE(size << ": "); for (auto it=begin;it!=end;it++) INFO_NONEWLINE(*it << " "); INFO("");
        return std::make_tuple(mean, variance, std::sqrt(variance/size));
        } 
};

/*
template <typename view_t>
void binning(const view_t& data, int max_bin_depth)
{
    view_t in(data);
    for (int i=0; i<max_bin_depth; i++) {
        size_t size = in.shape()[0];
        double mean, variance; std::tie(mean,variance) = calc_stats(in);
        MY_DEBUG(i <<": " << size << " " << mean << " " << variance);
        binning_adapter<view_t> b(in,i);
    }
}
*/

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_BINNING_HPP_

