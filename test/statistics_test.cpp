#include "common.hpp"
#include "triqs_extra.hpp"
#include "binning.hpp"
#include <triqs/mc_tools/random_generator.hpp>
#include <triqs/utility/tuple_tools.hpp>

using namespace fk;

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
    binning_adapter::binning<6>(a.begin(), a.end());
    }
catch (triqs::runtime_error const & e) { std::cerr  << "exception "<< e.what() << std::endl;}
   
    MY_DEBUG("--0 :" << binning_adapter::stats(a.begin(), a.end(), a.size()));
    MY_DEBUG("--0 :" << binning_adapter::stats(a.crbegin(), a.crend(), a.size()));
    MY_DEBUG(binning_adapter::binning<0>(a.begin(), a.end()));
    MY_DEBUG(binning_adapter::binning<0>(a.rbegin(), a.rend()));
    MY_DEBUG(binning_adapter::binning<1>(a.begin(), a.end()));
    MY_DEBUG(binning_adapter::binning<2>(a.begin(), a.end()));

    MY_DEBUG(binning(a.begin(), a.end(),2));
    MY_DEBUG(binning(a,2));

    auto bin_stats = accumulate_binning<3>::errors(a.begin(),a.end());
    MY_DEBUG(bin_stats.size());
    for (auto x : bin_stats) MY_DEBUG(x);
}
