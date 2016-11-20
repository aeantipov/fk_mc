#ifndef __FK_MC_BINNING_HPP_
#define __FK_MC_BINNING_HPP_

#include "common.hpp"

#include <cmath>
#include <ctgmath>
#include <boost/iterator/iterator_facade.hpp>

#include <boost/preprocessor.hpp>
#undef BOOST_PP_VARIADICS 
#define BOOST_PP_VARIADICS 1 
#include <iterator>

namespace fk {
namespace binning { 

#define BINNING_RANGE ((15)(14)(13)(12)(11)(10)(9)(8)(7)(6)(5)(4)(3)(2)(1)(0))

enum bin_m : size_t { _SIZE, _MEAN, _DISP, _SQERROR };
typedef std::tuple<size_t, double, double, double> bin_stats_t;
typedef std::vector<bin_stats_t> bin_data_t;


/** Binned iterator takes bin as an input parameter. If bin = 0 -> normal iterator, if bin>0 dereference returns averaged over (2^bin) values. */
template <typename iter_t, size_t D>
struct binned_iterator : 
    public boost::iterator_facade<binned_iterator<iter_t,D>, const double, boost::random_access_traversal_tag, double>
{
    typedef iter_t iterator_type;
    constexpr static size_t D_ = D;
    constexpr static size_t step = binned_iterator<iter_t,D-1>::step*2;
    iter_t p;
    binned_iterator(){};
    binned_iterator(iter_t start):p(start){};
    template <typename it2 = iter_t, typename std::enable_if<std::is_convertible<typename std::iterator_traits<it2>::iterator_category,std::random_access_iterator_tag>::value,bool>::type=0> 
    constexpr double dereference() const { return (*binned_iterator<iter_t,D-1>(p) + *binned_iterator<iter_t,D-1>(p+binned_iterator<iter_t,D-1>::step))/2.; };
    template <typename it2 = iter_t, typename std::enable_if<std::is_same<typename std::iterator_traits<it2>::iterator_category,std::forward_iterator_tag>::value,bool>::type=0> 
    double dereference() const { iter_t d(p); std::advance(d, binned_iterator<iter_t,D-1>::step); return (*binned_iterator<iter_t,D-1>(p) + *binned_iterator<iter_t,D-1>(d))/2.;};
    bool equal(binned_iterator const &rhs) const { return (p == rhs.p); };
    void increment() { std::advance(p,step); };
    void decrement() { std::advance(p,-int(step)); };
    void advance(int n) { std::advance(p,n*int(step)); };
    int distance_to(binned_iterator r){return std::distance(p,r.p)/step;};
};

template <typename iter_t>
struct binned_iterator<iter_t,0>: public boost::iterator_facade<binned_iterator<iter_t,0>, double, boost::random_access_traversal_tag, double>
{
    typedef iter_t iterator_type;
    constexpr static size_t D_ = 0;
    constexpr static size_t step = 1;
    iter_t p;
    binned_iterator(){};
    binned_iterator(iter_t start):p(start){};
    bool equal(binned_iterator const &rhs) const { return (p == rhs.p); };
    //template <typename std::enable_if<std::is_same<typename std::iterator_traits<iter_t>::iterator_category,std::random_access_iterator_tag>::value,bool>::type=0> 
    constexpr double dereference() const { return *p; };
    void increment() { p++; };
    void decrement() { p--; };
    void advance(int n) { std::advance(p,n); };
    int distance_to(binned_iterator r){return std::distance(p,r.p);};
};

template <size_t bin>
inline constexpr size_t bin_size(){return binned_iterator<std::vector<double>::iterator,bin>::size;};

inline size_t calc_bin_size(size_t size, size_t bin_step)
{
    if (bin_step > size) {
        FKMC_ERROR << "Can't bin with binning step(" << bin_step << ")> container size (" << size << ")";
        };
    int nsteps = size/bin_step;
    return nsteps;
}

template <typename bin_it, typename iter_t>
std::pair<bin_it,bin_it> find_bin_range(iter_t begin, iter_t end, bin_it)
{    
    int size = std::distance(begin, end);
    int step = bin_it::step;
    int nsteps = calc_bin_size(size,step);
    bin_it binned_end(begin); std::advance(binned_end,nsteps);
    return std::make_pair(bin_it(begin),binned_end);
}

//////

template <typename iter_t>
static bin_stats_t calc_stats(const iter_t& begin, const iter_t& end, size_t size = 0) 
{
    if(!size) size = std::distance(begin, end);
    double mean = std::accumulate(begin,end,0.0, std::plus<double>())/size; 
    double variance = std::accumulate(begin,end,0.0, [mean](double x, double y){return x + (y-mean)*(y-mean);})/(size-1);
    return std::make_tuple(size, mean, variance, std::sqrt(variance/size));
} 



template <size_t D, typename iter_t> static bin_stats_t bin(iter_t begin, iter_t end) 
{
    auto t1 = find_bin_range(begin,end,binned_iterator<iter_t, D>());
    return calc_stats(t1.first, t1.second);
};

template <size_t total_bins_left, size_t current_bin = 0>
struct binning_accumulator {
    template <typename iter_t>
    static bin_data_t accumulate_binning(iter_t begin, iter_t end) {
        auto stats = bin<current_bin>(begin,end); 
        auto next  = binning_accumulator<total_bins_left - 1, current_bin + 1>::accumulate_binning(begin,end);
        next.insert(next.begin(),stats);
        return next;
    };
};

template <size_t current_bin>
struct binning_accumulator<0,current_bin> {
    template <typename iter_t>
    static bin_data_t accumulate_binning(iter_t begin, iter_t end) {
        bin_data_t out;
        out.reserve(current_bin);
        auto stats = bin<current_bin>(begin,end);
        out.push_back(stats);
        return out;
    };
};


// free functions
template <typename iter_t>
static bin_stats_t bin(iter_t begin, iter_t end, int bin_depth) {
    #define MACRO(r, p) \
    if (BOOST_PP_SEQ_ELEM(0, p) == bin_depth) \
        return bin<BOOST_PP_SEQ_ELEM(0, p), iter_t>(begin,end);
    BOOST_PP_SEQ_FOR_EACH_PRODUCT(MACRO, BINNING_RANGE)
    #undef MACRO
    FKMC_ERROR << "bin_depth =" << bin_depth << "> compiled bin size";
    return std::make_tuple(0, std::nan(""),std::nan(""),std::nan(""));
    };

template <typename container_t>
static bin_stats_t bin(const container_t& in, int bin_depth) {
    return bin(in.begin(),in.end(),bin_depth);};

template <typename iter_t>
static bin_data_t accumulate_binning(iter_t begin, iter_t end, size_t bin_depth) {
    #define MACRO(r, p) \
    if (BOOST_PP_SEQ_ELEM(0, p) == bin_depth) \
        return binning_accumulator<BOOST_PP_SEQ_ELEM(0, p)>::accumulate_binning(begin,end);
    BOOST_PP_SEQ_FOR_EACH_PRODUCT(MACRO, BINNING_RANGE)
    #undef MACRO
    FKMC_ERROR << "bin_depth =" << bin_depth << "> compiled bin size";
    return {std::make_tuple(0, std::nan(""),std::nan(""),std::nan(""))};
}

template <typename container_t>
static bin_data_t accumulate_binning(const container_t& in, size_t bin_depth) {
    return accumulate_binning(in.begin(),in.end(),bin_depth);
}


inline std::vector<double> calc_cor_length(const bin_data_t& in)
{
    double sigma = std::get<_DISP>(in[0]);
    std::vector<double> out; out.reserve(in.size());
    for (size_t i=0; i<in.size(); ++i) {
        out[i] = 0.5*(std::pow(2.,i)*std::get<_DISP>(in[i])/sigma-1); 
    }
    return out;
}

#undef BINNING_RANGE

} // end of namespace binning
} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_BINNING_HPP_

