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
    std::cout << "Data A: [" << std::flush; for (auto x : a) std::cout << x << " " << std::flush; std::cout << "]" << std::endl;
    std::cout << "Data B: [" << std::flush; for (auto x : b) std::cout << x << " " << std::flush; std::cout << "]" << std::endl;

    typedef std::function<double(double)> xf_t;
    xf_t f1 = [](double x){return x;};

    /********************* part 1  **********************/
    INFO("T1: jackknife_adapter::jack test with trivial function");
    typedef std::vector<double>::iterator it_t;
    std::pair<it_t,it_t> orig_range = std::make_pair(a.begin(),a.end());
    //std::array<std::pair<it_t,it_t>,1> o = { orig_range };
    auto res = jackknife_adapter::jack<xf_t, it_t> (f1,{orig_range});
    bin_stats_t res_correct = std::make_tuple(17,0.484035,0.0872002,0.07162);
    MY_DEBUG(res << "==" << res_correct);
    if (!is_equal(res,res_correct,1e-5)) return EXIT_FAILURE;
    
    typedef binning::binned_iterator<it_t,0> bin0_it;
    auto res0 = jackknife_adapter::jack<xf_t,bin0_it>(f1,{find_bin_range(a.begin(),a.end(), bin0_it())});
    MY_DEBUG(res << "==" << res0);
    if (!is_equal(res0,res,1e-6)) return EXIT_FAILURE;

    typedef binning::binned_iterator<it_t,1> bin1_it;
    res = jackknife_adapter::jack<xf_t,bin1_it>(f1,{find_bin_range(a.begin(),a.end(), bin1_it())});
    res_correct = std::make_tuple(8,0.467231,0.0422793,0.0726974);
    MY_DEBUG(res << "==" << res_correct);
    if (!is_equal(res,res_correct,1e-5)) return EXIT_FAILURE;

    typedef binning::binned_iterator<it_t,2> bin2_it;
    res = jackknife_adapter::jack<xf_t,bin2_it>(f1,{find_bin_range(a.begin(),a.end(), bin2_it())});
    res_correct = std::make_tuple(4,0.467231,0.0269458,0.0820759);
    MY_DEBUG(res << "==" << res_correct);
    if (!is_equal(res,res_correct,1e-5)) return EXIT_FAILURE;

    /********************* part 2  **********************/
    INFO("T2: jackknife_adapter::jack test with nontrivial function");

    typedef std::function<double(double, double, double)> cf_t;
    cf_t f2 = [](double e, double e2, double de2){return e2 - de2 - e*e;}; 
    res = jackknife_adapter::jack<cf_t, it_t> (f2,{orig_range, std::make_pair(a2.begin(), a2.end()), std::make_pair(b.begin(),b.end())});
    res_correct = std::make_tuple(17,0.0297767,0.00815104,0.0218969);
    MY_DEBUG(res << "==" << res_correct);
    if (!is_equal(res,res_correct,1e-5)) return EXIT_FAILURE;

    res = jackknife_adapter::jack<cf_t, bin1_it> (f2,{
            find_bin_range(a.begin(),a.end(), bin1_it()),
            find_bin_range(a2.begin(),a2.end(), bin1_it()),
            find_bin_range(b.begin(),b.end(), bin1_it())
            });                
    res_correct = std::make_tuple(8,0.0309895,0.00214726,0.0163832);
    MY_DEBUG(res << "==" << res_correct);
    if (!is_equal(res,res_correct,1e-5)) return EXIT_FAILURE;

    res = jackknife_adapter::jack<cf_t, bin2_it> (f2,{
            find_bin_range(a.begin(),a.end(), bin2_it()),
            find_bin_range(a2.begin(),a2.end(), bin2_it()),
            find_bin_range(b.begin(),b.end(), bin2_it())
            });                
    res_correct = std::make_tuple(4,0.0324411,0.00123011,0.0175365);
    MY_DEBUG(res << "==" << res_correct);
    if (!is_equal(res,res_correct,1e-5)) return EXIT_FAILURE;

    /********************* part 3  **********************/
    INFO("T3: accumulate_jackknife"); 
    auto stats = accumulate_jackknife(f1,std::vector<std::pair<it_t,it_t>>({orig_range}),3);
    for (auto x:stats){MY_DEBUG(x);}; 
    res = stats[3];
    res_correct = std::make_tuple(2,0.467231,0.00162846,0.0285348);
    MY_DEBUG(res << "==" << res_correct);
    if (!is_equal(res,res_correct,1e-5)) return EXIT_FAILURE;

    auto stats2 = accumulate_jackknife(f1,std::vector<std::vector<double>>({a}),3);
    for (auto x:stats){MY_DEBUG(x);}; 
    for (size_t i=0; i<stats.size(); i++) if (!is_equal(stats[i],stats2[i],1e-5)) return EXIT_FAILURE;

    /********************* part 4  **********************/

    INFO("T3: accumulate_jackknife with a function with a lot of arguments"); 
    constexpr size_t nargs = 12; // need f(x_1, x_2, ... x_12);
    size_t nmeasures = 17;
    std::vector<std::vector<double>> data(nargs);
    for (size_t i=0; i<nargs; i++) { data[i].resize(nmeasures); for (size_t j=0; j<nmeasures; j++) data[i][j]=RNG(); };
    for (size_t j=0; j<nmeasures; j++) {
        for (size_t i=0; i<nargs; i++)  
            std::cout << data[i][j] <<" " << std::flush; std::cout << std::endl; 
        };

    std::function<double(double,double,double,double,double,double,double,double,double,double,double,double)> f4_1 = 
               [](double x1, double x2, double x3, double x4, double x5, double x6, double x7, double x8, double x9, double x10, double x11, double x12){return x1+x2;};
    stats = accumulate_jackknife(f4_1,data,2);
    for (auto x:stats){MY_DEBUG(x);}; 
    
    bool success = true;
    if (!success) return EXIT_FAILURE;
    return EXIT_SUCCESS;
}
