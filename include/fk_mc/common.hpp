#ifndef __FK_MC_COMMON_H
#define __FK_MC_COMMON_H

#include "definitions.hpp"

//#include <triqs/arrays.hpp>
//#include <triqs/arrays/algorithms.hpp>
//#include <triqs/parameters/parameters.hpp>
#include <complex>

#include <tuple>
#include <array>
#include <boost/math/special_functions/pow.hpp>
#include <numeric>
#include <random>
#include <alps/params.hpp>

namespace fk {

typedef std::mt19937 random_generator;

//namespace tqa = triqs::arrays;

//using triqs::arrays::prod;
//using triqs::arrays::sum;

// FKDEBUG messages with custom verbosity.
// Adapted from http://efesx.com/2010/08/31/overloading-macros/
#ifndef NDEBUG
#define MSG_PREFIX            __FILE__ << ":" << __LINE__ << ": "
#define FKDEBUG3(MSG,VERBOSITY,VERB_LEVEL)            if (VERBOSITY >= VERB_LEVEL) std::cerr << MSG_PREFIX << MSG << std::endl;
#else
#define FKDEBUG3(MSG,VERBOSITY,VERB_LEVEL)            ;
#endif

#define FKDEBUG1(MSG) FKDEBUG3(MSG,3,3)
#define FKDEBUG2(MSG,VERBOSITY) FKDEBUG3(MSG,VERBOSITY,3)

#define VA_NUM_ARGS_IMPL(_1,_2,_3,_4,_5,N,...) N
#define VA_NUM_ARGS(...) VA_NUM_ARGS_IMPL(__VA_ARGS__, 5,4,3,2,1)
#define macro_dispatcher(func, ...) \
            macro_dispatcher_(func, VA_NUM_ARGS(__VA_ARGS__))
#define macro_dispatcher_(func, nargs) \
            macro_dispatcher__(func, nargs)
#define macro_dispatcher__(func, nargs) \
            func ## nargs

#define FKDEBUG(...) macro_dispatcher(FKDEBUG, __VA_ARGS__)(__VA_ARGS__)


/*
#define INFO(MSG)             std::cout << std::boolalpha << MSG << std::endl;
#define INFO2(MSG)            std::cout << "    " << std::boolalpha << MSG << std::endl;
#define INFO_NONEWLINE(MSG)   std::cout << MSG << std::flush;
*/
//inline std::ostream& err_and_throw(std::string x, std::ostream& y) { y << x; throw std::logic_error(x); } 
struct err_and_throw { 
    std::ostream& str_;
    template <typename T>
    friend err_and_throw& operator<< (err_and_throw&& y, T x) { y.str_ << x; std::stringstream s1; s1 << x; throw std::logic_error(s1.str()); return y; } 
    template <typename T>
    friend err_and_throw& operator<< (err_and_throw& y, T x) { y.str_ << x; std::stringstream s1; s1 << x; throw std::logic_error(s1.str()); return y; } 
    err_and_throw(std::ostream &x):str_(x) {}
    //static err_and_throw &get_instance() { static err_and_throw x; return x; };
    };
#define TRIQS_RUNTIME_ERROR err_and_throw(std::cerr)

//typedef triqs::utility::parameters parameters_t;
typedef alps::params parameters_t;
typedef std::complex<double> complex_t;
static const complex_t I (0.0,1.0); 
static const double PI = atan(1)*4.;

template <size_t D>
inline std::ostream& operator<< (std::ostream& in, std::array<int, D> arr){ 
    in << "{";  std::copy(arr.begin(), arr.end()-1, std::ostream_iterator<int>(in, ",")); in << *(arr.end()-1) << "}";
    return in;
};

/** A tool to split a tuple from http://stackoverflow.com/questions/10626856/how-to-split-a-tuple. */
template <typename Target, typename Tuple, int N, bool end >
struct __split_tuple_struct
{
    template < typename ... Args >
    static Target create(Tuple const& t, Args && ... args)
    {
        return __split_tuple_struct<Target,Tuple, N+1, std::tuple_size<Tuple>::value == N+1>::create(t, std::forward<Args>(args)..., std::get<N>(t));
    }
};

template < typename Target, typename Tuple, int N >
struct __split_tuple_struct<Target,Tuple,N,true>
{
    template < typename ... Args >
    static Target create(Tuple const& t, Args && ... args) { return Target(std::forward<Args>(args)...); }
};

template < typename Head, typename ... Tail >
inline std::tuple<Tail...> tuple_tail(std::tuple<Head,Tail...> const& tpl)
{
    return __split_tuple_struct<std::tuple<Tail...>, std::tuple<Head,Tail...>, 1, std::tuple_size<std::tuple<Head,Tail...>>::value == 1>::create(tpl);
}

// function_traits from http://stackoverflow.com/questions/9044866/how-to-get-the-number-of-arguments-of-stdfunction

template<typename T> struct function_traits;
template<typename R, typename ...Args> 
struct function_traits<std::function<R(Args...)>>
{
    static constexpr size_t nargs = sizeof...(Args);
    typedef R result_type;
    template <size_t i>
    struct arg { typedef typename std::tuple_element<i, std::tuple<Args...>>::type type; };
};

}; // end of namespace FK

#endif // endif :: #ifndef __FK_MC_COMMON_H


