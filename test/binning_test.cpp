#include "common.hpp"
#include "binning.hpp"
#include <gftools/tuple_tools.hpp>

using namespace fk;
using namespace binning;

template <typename T1> struct float_comparator {
static bool is_equal(T1 a, T1 b, double tol){return (std::abs(double(a)-double(b))<tol);};
};

template <typename Arg1> struct float_comparator<std::tuple<Arg1>> {
static bool is_equal(std::tuple<Arg1> a, std::tuple<Arg1> b, double tol){return float_comparator<Arg1>::is_equal(std::get<0>(a), std::get<0>(b), tol);};
};


template <typename Arg1, typename ... Args> struct float_comparator<std::tuple<Arg1, Args...>> {
static bool is_equal (std::tuple<Arg1, Args...> a, std::tuple<Arg1, Args...> b, double tol)
{ 
return (float_comparator<Arg1>::is_equal(std::get<0>(a),std::get<0>(b),tol) && float_comparator<std::tuple<Args...>>::is_equal(tuple_tail(a),tuple_tail(b), tol));};
};

template <typename T1> bool is_equal(T1 a, T1 b, double tol = std::numeric_limits<double>::epsilon()){return float_comparator<T1>::is_equal(a,b,tol);};

using namespace gftools::tuple_tools;

int main()
{

    std::vector<double> a = {0.0711992, 0.344949, 0.940913, 0.166604, 0.811305, 0.617859, 
0.462844, 0.550449, 0.28126, 0.0560575, 0.0673806, 0.710085, 
0.459742, 0.977218, 0.500193, 0.45763, 0.752903};

    binned_iterator<std::vector<double>::iterator,0> it0(a.begin());
    FKDEBUG(it0.step); FKDEBUG(*it0); it0++; FKDEBUG(*it0);

    binned_iterator<std::vector<double>::iterator,1> it1(a.begin());
    FKDEBUG(it1.step); FKDEBUG(*it1); it1++; FKDEBUG(*it1);

    binned_iterator<std::vector<double>::iterator,2> it2(a.begin());
    FKDEBUG(it2.step); FKDEBUG(*it2); it2++; FKDEBUG(*it2);

try{ 
    bin<6>(a.begin(), a.end());
    }
catch (std::exception const & e) { std::cerr  << "exception "<< e.what() << std::endl;}
   
    FKDEBUG("--0 :" << calc_stats(a.begin(), a.end()));
    FKDEBUG("--0 :" << calc_stats(a.crbegin(), a.crend()));
    FKDEBUG(bin<0>(a.begin(), a.end()));
    FKDEBUG(bin<0>(a.rbegin(), a.rend()));
    FKDEBUG(bin<1>(a.begin(), a.end()));
    FKDEBUG(bin<2>(a.begin(), a.end()));

    FKDEBUG(bin(a.begin(), a.end(),2));
    FKDEBUG(bin(a,2));

    auto bin_stats = binning_accumulator<3>::accumulate_binning(a.begin(),a.end());
    auto cor_lens = calc_cor_length(bin_stats); 
    std::vector<std::tuple<size_t,double,double,double>> correct_v = // taken from Mathematica
    { std::make_tuple(17, 0.484035, 0.0872001, 0.07162), std::make_tuple(8, 0.467231, 0.0422793, 0.0726974), 
      std::make_tuple(4, 0.467231, 0.0269458, 0.0820759), std::make_tuple(2, 0.467231,0.00162846,  0.0285348)};
    std::vector<double> correct_cor_lens({0, -0.0151463, 0.118023, -0.4253});
    bool success = true;
    for (size_t i=0; i<bin_stats.size(); i++) {
        INFO(i <<" : " << print_tuple(bin_stats[i]) << "," << cor_lens[i] << "==" << print_tuple(correct_v[i]) << "," << correct_cor_lens[i] << 
            " = " << is_equal(bin_stats[i],correct_v[i],1e-5) << " " << is_equal(cor_lens[i],correct_cor_lens[i],1e-5));
        success = success && is_equal(bin_stats[i],correct_v[i],1e-5) && is_equal(cor_lens[i],correct_cor_lens[i],1e-5); 
    };
    if (!success) return EXIT_FAILURE;

    auto bin_stats2 = accumulate_binning(a.begin(),a.end(), 3); // call non-templated
    for (size_t i=0; i<bin_stats2.size(); i++) {
        INFO(i <<" : " << print_tuple(bin_stats2[i]) << "," << cor_lens[i] << "==" << print_tuple(correct_v[i]) << "," << correct_cor_lens[i] << 
            " = " << is_equal(bin_stats2[i],correct_v[i],1e-5) << " " << is_equal(cor_lens[i],correct_cor_lens[i],1e-5));
        success = success && is_equal(bin_stats2[i],correct_v[i],1e-5) && is_equal(cor_lens[i],correct_cor_lens[i],1e-5); 
    };
    if (!success) return EXIT_FAILURE;


    std::array<double, 17> b; std::copy(a.begin(),a.end(),b.begin());
    auto bin_stats_tqa = binning_accumulator<3>::accumulate_binning(b.begin(),b.end());
    for (size_t i=0; i<bin_stats_tqa.size(); i++) {
        INFO(i <<" : " << print_tuple(bin_stats_tqa[i]) << "," << cor_lens[i] << "==" << print_tuple(correct_v[i]) << "," << correct_cor_lens[i] << 
            " = " << is_equal(bin_stats_tqa[i],correct_v[i],1e-5) << " " << is_equal(cor_lens[i],correct_cor_lens[i],1e-5));
        success = success && is_equal(bin_stats_tqa[i],correct_v[i],1e-5) && is_equal(cor_lens[i],correct_cor_lens[i],1e-5); 
    };
    if (!success) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
