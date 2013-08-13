#include "common.hpp"
#include "triqs_extra.hpp"
#include "jackknife.hpp"
#include <triqs/mc_tools/random_generator.hpp>
#include <triqs/utility/tuple_tools.hpp>

using namespace fk;
using namespace jackknife;

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
    std::vector<double> a = {0.0711992, 0.344949, 0.940913, 0.166604, 0.811305, 0.617859, 
0.462844, 0.550449, 0.28126, 0.0560575, 0.0673806, 0.710085, 
0.459742, 0.977218, 0.500193, 0.45763, 0.752903};
    std::vector<double> a2 (a.size());
    std::transform(a.begin(), a.end(), a2.begin(), [](double x){ return x*x; });
    std::vector<double> b = {0.0203252 , 0.0491541, 0.0942537, 0.0815104, 0.0569245, 0.0459406, 
  0.0963107, 0.0170473, 0.0730589, 0.012731, 0.0571613, 0.0889872, 
  0.0771268, 0.0857561, 0.0130207, 0.0378117, 0.0690792};
    //auto view1 = fk::make_weak_view(a);
    tqa::vector<double> tqa_a(a.size());
    std::copy(a.begin(),a.end(),tqa_a.begin());
    tqa::vector<double> tqa_b(b.size());
    std::copy(b.begin(),b.end(),tqa_b.begin());
    MY_DEBUG(tqa_a);
    MY_DEBUG(tqa_b);

    typedef std::function<double(double)> xf_t;
    xf_t f1 = [](double x){return x;};

    typedef std::vector<double>::iterator it_t;
    std::pair<it_t,it_t> orig_range = std::make_pair(a.begin(),a.end());
    //std::array<std::pair<it_t,it_t>,1> o = { orig_range };
    auto res = jackknife_adapter::jack<xf_t, it_t, 1> (f1,{orig_range});
    MY_DEBUG(res);
    
    typedef binning::binned_iterator<it_t,0> bin0_it;
    res = jackknife_adapter::jack<xf_t,bin0_it,1>(f1,{find_bin_range(a.begin(),a.end(), bin0_it())});
    MY_DEBUG(res);

    typedef binning::binned_iterator<it_t,1> bin1_it;
    res = jackknife_adapter::jack<xf_t,bin1_it,1>(f1,{find_bin_range(a.begin(),a.end(), bin1_it())});
    MY_DEBUG(res);

    typedef binning::binned_iterator<it_t,2> bin2_it;
    res = jackknife_adapter::jack<xf_t,bin2_it,1>(f1,{find_bin_range(a.begin(),a.end(), bin2_it())});
    MY_DEBUG(res);


    typedef std::function<double(double, double, double)> cf_t;
    cf_t f2 = [](double e, double e2, double de2){return e2 - de2 - e*e;}; 
    res = jackknife_adapter::jack<cf_t, it_t, 3> (f2,{orig_range, std::make_pair(a2.begin(), a2.end()), std::make_pair(b.begin(),b.end())});
    MY_DEBUG(res);
    res = jackknife_adapter::jack<cf_t, bin1_it, 3> (f2,{
            find_bin_range(a.begin(),a.end(), bin1_it()),
            find_bin_range(a2.begin(),a2.end(), bin1_it()),
            find_bin_range(b.begin(),b.end(), bin1_it())
            });                
    MY_DEBUG(res);
    res = jackknife_adapter::jack<cf_t, bin2_it, 3> (f2,{
            find_bin_range(a.begin(),a.end(), bin2_it()),
            find_bin_range(a2.begin(),a2.end(), bin2_it()),
            find_bin_range(b.begin(),b.end(), bin2_it())
            });                
    MY_DEBUG(res);

    auto stats = accumulate_jackknife(f1,std::array<std::pair<it_t,it_t>, 1>({orig_range}),3);
    for (auto x:stats){MY_DEBUG(x);}; 
    stats = accumulate_jackknife(f1,std::array<std::vector<double>, 1>({a}),3);
    for (auto x:stats){MY_DEBUG(x);}; 

    std::vector<double> d(30000);
    for (auto &x : d) x = RNG();
    auto d_range = std::make_pair(d.begin(), d.end());
    //accumulate_jackknife(f1,std::array<std::pair<it_t,it_t>, 1>({std::make_pair(d.begin(), d.end())}),10);
    res = jackknife_adapter::jack<xf_t, it_t, 1> (f1,{d_range});
    MY_DEBUG(res);

    
    bool success = true;
    if (!success) return EXIT_FAILURE;
    return EXIT_SUCCESS;
}
