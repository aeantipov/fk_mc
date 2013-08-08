#include "common.hpp"
#include "triqs_extra.hpp"
#include "binning.hpp"
#include <triqs/mc_tools/random_generator.hpp>
#include <triqs/utility/tuple_tools.hpp>

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

int main()
{
    triqs::mc_tools::random_generator RNG("mt19937", 23432);
    std::vector<double> a(17,1.5);
    for (auto &b: a) { b= RNG() ;std::cout << b <<" " << std::flush; }; INFO("");
    //auto view1 = fk::make_weak_view(a);
    tqa::vector<double> arr1(a.size());
    std::copy(a.begin(),a.end(),arr1.begin());
    
    INFO(arr1());

    binned_iterator<std::vector<double>::iterator,0> it0(a.begin());
    MY_DEBUG(it0.step); MY_DEBUG(*it0); it0++; MY_DEBUG(*it0);

    binned_iterator<std::vector<double>::iterator,1> it1(a.begin());
    MY_DEBUG(it1.step); MY_DEBUG(*it1); it1++; MY_DEBUG(*it1);

    binned_iterator<std::vector<double>::iterator,2> it2(a.begin());
    MY_DEBUG(it2.step); MY_DEBUG(*it2); it2++; MY_DEBUG(*it2);

try{ 
    binning_adapter::bin<6>(a.begin(), a.end());
    }
catch (triqs::runtime_error const & e) { std::cerr  << "exception "<< e.what() << std::endl;}
   
    MY_DEBUG("--0 :" << binning_adapter::calc_stats(a.begin(), a.end()));
    MY_DEBUG("--0 :" << binning_adapter::calc_stats(a.crbegin(), a.crend()));
    MY_DEBUG(binning_adapter::bin<0>(a.begin(), a.end()));
    MY_DEBUG(binning_adapter::bin<0>(a.rbegin(), a.rend()));
    MY_DEBUG(binning_adapter::bin<1>(a.begin(), a.end()));
    MY_DEBUG(binning_adapter::bin<2>(a.begin(), a.end()));

    MY_DEBUG(bin(a.begin(), a.end(),2));
    MY_DEBUG(bin(a,2));

    auto bin_stats = binning_accumulator<3>::accumulate_binning(a.begin(),a.end());
    auto cor_lens = calc_cor_length(bin_stats); 
    std::vector<std::tuple<size_t,double,double,double>> correct_v = // taken from Mathematica
    { std::make_tuple(17, 0.484035, 0.0872001, 0.07162), std::make_tuple(8, 0.467231, 0.0422793, 0.0726974), 
      std::make_tuple(4, 0.467231, 0.0269458, 0.0820759), std::make_tuple(2, 0.467231,0.00162846,  0.0285348)};
    std::vector<double> correct_cor_lens({0, -0.0151463, 0.118023, -0.4253});
    bool success = true;
    for (size_t i=0; i<bin_stats.size(); i++) {
        INFO(i <<" : " << bin_stats[i] << "," << cor_lens[i] << "==" << correct_v[i] << "," << correct_cor_lens[i] << 
            " = " << is_equal(bin_stats[i],correct_v[i],1e-5) << " " << is_equal(cor_lens[i],correct_cor_lens[i],1e-5));
        success = success && is_equal(bin_stats[i],correct_v[i],1e-5) && is_equal(cor_lens[i],correct_cor_lens[i],1e-5); 
    };
    if (!success) return EXIT_FAILURE;

    auto bin_stats2 = accumulate_binning(a.begin(),a.end(), 3); // call non-templated
    for (size_t i=0; i<bin_stats2.size(); i++) {
        INFO(i <<" : " << bin_stats2[i] << "," << cor_lens[i] << "==" << correct_v[i] << "," << correct_cor_lens[i] << 
            " = " << is_equal(bin_stats2[i],correct_v[i],1e-5) << " " << is_equal(cor_lens[i],correct_cor_lens[i],1e-5));
        success = success && is_equal(bin_stats2[i],correct_v[i],1e-5) && is_equal(cor_lens[i],correct_cor_lens[i],1e-5); 
    };
    if (!success) return EXIT_FAILURE;


    tqa::array<double,1> b(17); std::copy(a.begin(),a.end(),b.begin());
    auto bin_stats_tqa = binning_accumulator<3>::accumulate_binning(b.begin(),b.end());
    for (size_t i=0; i<bin_stats_tqa.size(); i++) {
        INFO(i <<" : " << bin_stats_tqa[i] << "," << cor_lens[i] << "==" << correct_v[i] << "," << correct_cor_lens[i] << 
            " = " << is_equal(bin_stats_tqa[i],correct_v[i],1e-5) << " " << is_equal(cor_lens[i],correct_cor_lens[i],1e-5));
        success = success && is_equal(bin_stats_tqa[i],correct_v[i],1e-5) && is_equal(cor_lens[i],correct_cor_lens[i],1e-5); 
    };
    if (!success) return EXIT_FAILURE;

    return EXIT_SUCCESS;
}
