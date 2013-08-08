#ifndef __FK_MC_BINNING_HPP_
#define __FK_MC_BINNING_HPP_

#include "common.hpp"
#include "triqs_extra.hpp"

#include <cmath>
#include <ctgmath>
#include <boost/iterator/iterator_facade.hpp>
#include <triqs/utility/exceptions.hpp>
#include <boost/preprocessor.hpp>
#include <iterator>

namespace fk {

#define BINNING_RANGE ((15)(14)(13)(12)(11)(10)(9)(8)(7)(6)(5)(4)(3)(2)(1)(0))

/** Binned iterator takes bin as an input parameter. If bin = 0 -> normal iterator, if bin>0 dereference returns averaged over (2^bin) values. */
template <typename iter_t, size_t D>
struct binned_iterator : 
    public boost::iterator_facade<binned_iterator<iter_t,D>, const double, boost::random_access_traversal_tag, double>
{
    constexpr static size_t step = binned_iterator<iter_t,D-1>::step*2;
    iter_t p;
    binned_iterator(iter_t start):p(start){};
    template <typename it2 = iter_t, typename std::enable_if<std::is_convertible<typename std::iterator_traits<it2>::iterator_category,std::random_access_iterator_tag>::value,bool>::type=0> 
    constexpr double dereference() const { return (*binned_iterator<iter_t,D-1>(p) + *binned_iterator<iter_t,D-1>(p+binned_iterator<iter_t,D-1>::step))/2.; };
    template <typename it2 = iter_t, typename std::enable_if<std::is_same<typename std::iterator_traits<it2>::iterator_category,std::forward_iterator_tag>::value,bool>::type=0> 
    double dereference() const { iter_t d(p); std::advance(d, binned_iterator<iter_t,D-1>::step); return (*binned_iterator<iter_t,D-1>(p) + *binned_iterator<iter_t,D-1>(d))/2.;};
    bool equal(binned_iterator const &rhs) const { return (p == rhs.p); };
    void increment() { std::advance(p,step); };
};

template <typename iter_t>
struct binned_iterator<iter_t,0>: public boost::iterator_facade<binned_iterator<iter_t,0>, double, boost::random_access_traversal_tag, double>
{
    constexpr static size_t step = 1;
    iter_t p;
    binned_iterator(iter_t start):p(start){};
    bool equal(binned_iterator const &rhs) const { return (p == rhs.p); };
    //template <typename std::enable_if<std::is_same<typename std::iterator_traits<iter_t>::iterator_category,std::random_access_iterator_tag>::value,bool>::type=0> 
    constexpr double dereference() const { return *p; };
    void increment() { p++; };
};

struct binning_adapter {
    template <size_t bin, typename iter_t> static std::tuple<size_t, double, double, double> binning(iter_t begin, iter_t end) {
        int size = std::distance(begin, end);
        int step = binned_iterator<iter_t,bin>::step;
        if (step > size) TRIQS_RUNTIME_ERROR << "Can't bin with depth = " << bin << ", binning step(" << step << ")> container size (" << size << ")";
        int nsteps = size/step;
        binned_iterator<iter_t,bin> binned_begin(begin), binned_end(begin); std::advance(binned_end,nsteps);
        return stats(binned_begin, binned_end);
        };

    template <typename iter_t>
    static std::tuple<size_t, double, double, double> stats(const iter_t& begin, const iter_t& end) {
        int size = std::distance(begin, end);
        double mean = std::accumulate(begin,end,0.0, std::plus<double>())/size; 
        double variance = std::accumulate(begin,end,0.0, [mean](double x, double y){return x + (y-mean)*(y-mean);})/(size-1);
        //INFO_NONEWLINE(size << ": "); for (auto it=begin;it!=end;it++) INFO_NONEWLINE(*it << " "); INFO("");
        return std::make_tuple(size, mean, variance, std::sqrt(variance/size));
        } 
    //template <typename iter_t>
    //static std::
};


template <size_t total_bins_left, size_t current_bin = 0>
struct binning_accumulator {
    template <typename iter_t>
    static std::vector<std::tuple<size_t,double,double,double,double>> accumulate_binning(iter_t begin, iter_t end, double sigma = 0) {
        size_t size; double mean, disp, sqerror;
        auto stats = binning_adapter::binning<current_bin>(begin,end); std::tie(size, mean, disp, sqerror) = stats; 
        if (!current_bin) {sigma = disp;};
        auto next = binning_accumulator<total_bins_left - 1, current_bin + 1>::accumulate_binning(begin,end,sigma);
        next.insert(next.begin(),std::make_tuple(size, mean,disp,sqerror,0.5*(boost::math::pow<current_bin>(2.)*disp/sigma-1)));
        return next;
    };
};

template <size_t current_bin>
struct binning_accumulator<0,current_bin> {
    template <typename iter_t>
    static std::vector<std::tuple<size_t,double,double,double,double>> accumulate_binning(iter_t begin, iter_t end, double sigma) {
        std::vector<std::tuple<size_t,double,double,double,double>> out;
        out.reserve(current_bin);
        size_t size; double mean, disp, sqerror;
        auto stats = binning_adapter::binning<current_bin>(begin,end); std::tie(size, mean, disp, sqerror) = stats; 
        out.push_back(std::make_tuple(size,mean,disp,sqerror,0.5*(boost::math::pow<current_bin>(2.)*disp/sigma-1)));
        return out;
    };
};


// free functions
template <typename iter_t>
static std::tuple<size_t,double,double,double> binning(iter_t begin, iter_t end, int bin_depth) {
    #define MACRO(r, p) \
    if (BOOST_PP_SEQ_ELEM(0, p) == bin_depth) \
        return binning_adapter::binning<BOOST_PP_SEQ_ELEM(0, p), iter_t>(begin,end);
    BOOST_PP_SEQ_FOR_EACH_PRODUCT(MACRO, BINNING_RANGE)
    #undef MACRO
    TRIQS_RUNTIME_ERROR << "bin_depth =" << bin_depth << "> compiled bin size";
    return std::make_tuple(0, std::nan(""),std::nan(""),std::nan(""));
    };

template <typename container_t>
static std::tuple<size_t,double,double,double> binning(const container_t& in, int bin_depth) {
    return binning(in.begin(),in.end(),bin_depth);};

template <typename iter_t>
static std::vector<std::tuple<size_t,double,double,double,double>> accumulate_binning(iter_t begin, iter_t end, size_t bin_depth) {
    double sigma;
    #define MACRO(r, p) \
    if (BOOST_PP_SEQ_ELEM(0, p) == bin_depth) \
        return binning_accumulator<BOOST_PP_SEQ_ELEM(0, p)>::accumulate_binning(begin,end,sigma);
    BOOST_PP_SEQ_FOR_EACH_PRODUCT(MACRO, BINNING_RANGE)
    #undef MACRO
    TRIQS_RUNTIME_ERROR << "bin_depth =" << bin_depth << "> compiled bin size";
    return {std::make_tuple(0, std::nan(""),std::nan(""),std::nan(""),std::nan(""))};
}

#undef BINNING_RANGE

} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_BINNING_HPP_

