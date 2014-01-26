#ifndef __FK_MC_JACKKNIFE_HPP_
#define __FK_MC_JACKKNIFE_HPP_

#include "binning.hpp"

namespace fk {
namespace jackknife { 

// defines maximum binning depth
#define BINNING_RANGE ((15)(14)(13)(12)(11)(10)(9)(8)(7)(6)(5)(4)(3)(2)(1)(0))

using namespace binning;

/** A jackknife_adapter does statistics analysis of an indirect measure of the data, that is
 * of a C=F(x_1,x_2,...,x_n), where x variables are measured directly and F is a given function
 */
struct jackknife_adapter
{

    /** A tool to convert a function of arguments to a function of vector of arguments, i.e.
     * f(x_1, x_2, ...) -> f([x_1, x_2, ...])
     */
    template <typename R, typename Arg1, typename ... Args>
    static std::function<R(std::vector<Arg1>)> _vectorize_f(const std::function<R(Arg1, Args...)>& F_in) {
        constexpr size_t L = sizeof...(Args)+1;
        std::function<R(std::vector<Arg1>)> out_;
        out_ = [&](const std::vector<Arg1>& in)->R { 
            if (in.size()!=L) TRIQS_RUNTIME_ERROR << "Argument size mismatch for jackknife";
            std::array<Arg1, L> p; std::copy(in.begin(), in.end(), p.begin());
            return triqs::tuple::apply(F_in,p);
            };
        MY_DEBUG("Converting f of " << L << " args to f(std::vector)");
        return out_;
        };

    /** An overload for  _vectorize_f to pass a function of vector untouched. */
    template <typename R, typename Arg1>
    static std::function<R(std::vector<Arg1>)> _vectorize_f(const std::function<R(std::vector<Arg1>)>& F_in) {MY_DEBUG("Passing f"); return F_in;};

    /** Run a jackknife analysis for a given data.
     * \param[in] F - a function that gives the required expression, i.e. F(x_1, x_2 ... )
     * \param[in] data - a vector of ranges of variables [x_1, x_2]
     */
    template <typename Functor, typename iter_t>
    static bin_stats_t jack(Functor F, const std::vector<std::pair<iter_t,iter_t>> &data )
    {
        typedef std::pair<iter_t,iter_t> iter_pair_t;

        size_t size = std::distance(data[0].first,data[0].second);
        size_t L = data.size();
        std::vector<double> x_means(L);
        for (size_t i=0; i<L; ++i) x_means[i] = std::get<_MEAN>(calc_stats(data[i].first, data[i].second));
        auto F2 = _vectorize_f(F);
        double U_0 = F2(x_means);//triqs::tuple::apply(F,x_means); 

        std::vector<double> U(size);
        std::vector<iter_t> iters(L); 
        for (size_t i=0; i<L; ++i) { iters[i]=data[i].first; };

        std::vector<double> xy_means(L); // averages without i value
        for (size_t j=0; j<size; ++j) { 
        xy_means.assign(L,0.0);
            for (size_t i=0; i<L; ++i) { 
                xy_means[i] = (double(size)*x_means[i]-*iters[i])/(double(size-1));
                iters[i]++;
                }
            double U_j = F2(xy_means);//triqs::tuple::apply(F,xy_means);
            U[j] = U_j;
            };
        auto U_stats = calc_stats(U.begin(),U.end());
        double U_bar = std::get<_MEAN>(U_stats);
        double U_average = U_0 - (size-1)*(U_bar - U_0);
        double dU = (size-1)*std::get<_SQERROR>(U_stats);
        double var_U = dU*dU*size; 
        //MY_DEBUG(U_0 << " " << U_average << "+/-" << dU);
        return std::make_tuple(size, U_average, var_U, dU);
    }
};

/** A structure to accumulate jackknife stats for different bins of provided data. */
template <class Functor, size_t total_bins_left, size_t current_bin = 0>
struct jackknife_accumulator {
    /** Accumulates jackknife data for given data and a function over this data. Returns a vector of binned data. */
    template <typename iter_t>
    static bin_data_t accumulate_jackknife(Functor f, const std::vector<std::pair<iter_t,iter_t>> &data) {
        size_t L = data.size();
        typedef binning::binned_iterator<iter_t,current_bin> bin_it;
        std::vector<std::pair<bin_it, bin_it>> range(L);
        for (size_t i=0; i<L; i++) range[i]=find_bin_range(data[i].first, data[i].second, bin_it());
        auto stats = jackknife_adapter::jack<Functor,bin_it>(f,range);
        auto next  = jackknife_accumulator<Functor, total_bins_left - 1, current_bin + 1>::accumulate_jackknife(f,data);
        next.insert(next.begin(),stats);
        return next;
    };
};

template <class Functor, size_t current_bin>
struct jackknife_accumulator<Functor,0,current_bin> {
    template <typename iter_t>
    static bin_data_t accumulate_jackknife(Functor f, const std::vector<std::pair<iter_t,iter_t>> &data) {
        size_t L = data.size();
        bin_data_t out;
        out.reserve(current_bin);
        typedef binning::binned_iterator<iter_t,current_bin> bin_it;
        std::vector<std::pair<bin_it, bin_it>> range(L);
        for (size_t i=0; i<L; i++) range[i]=find_bin_range(data[i].first, data[i].second, bin_it());
        auto stats = jackknife_adapter::jack<Functor,bin_it>(f,range);
        out.push_back(stats);
        return out;
    };
};

//////// free functions ///////////

/** Accumulate jackknife data for a given data (vector of ranges), function F and required maximum bin_depth. */
template <typename Functor, typename iter_t>
bin_data_t accumulate_jackknife(Functor F, const std::vector<std::pair<iter_t,iter_t>>& in, size_t bin_depth) {
    #define MACRO(r, p) \
    if (BOOST_PP_SEQ_ELEM(0, p) == bin_depth) \
        { return jackknife_accumulator<Functor,BOOST_PP_SEQ_ELEM(0, p)>::template accumulate_jackknife<iter_t>(F,in); };
    BOOST_PP_SEQ_FOR_EACH_PRODUCT(MACRO, BINNING_RANGE)
    #undef MACRO
    TRIQS_RUNTIME_ERROR << "bin_depth =" << bin_depth << "> compiled bin size";
    return {std::make_tuple(0, std::nan(""),std::nan(""),std::nan(""))};
}

/** Accumulate jackknife data for a given data (vector of containers), function F and required maximum bin_depth. */
template <typename Functor, typename Std_container_t>
bin_data_t accumulate_jackknife(Functor &&F, const std::vector<Std_container_t>& in, size_t bin_depth) {
    typedef typename std::add_const<Std_container_t>::type::const_iterator iter_t;
    size_t L = in.size();
    std::vector<std::pair<iter_t,iter_t>> ranges(L);
    for (size_t i=0; i<L; ++i) ranges[i]=std::make_pair(in[i].begin(),in[i].end());
    return accumulate_jackknife(F,ranges,bin_depth);
}

} // end of namespace jackknife
} // end of namespace fk

#endif // endif :: #ifndef __FK_MC_JACKKNIFE_HPP_

