#ifndef __FK_MC_COMMON_H
#define __FK_MC_COMMON_H

#include <triqs/arrays.hpp>
#include <triqs/arrays/algorithms.hpp>
#include <triqs/parameters/parameters.hpp>

#include <tuple>
#include <array>
#include <boost/math/special_functions/pow.hpp>
#include <numeric>
#include<Eigen/Core>

namespace fk {

namespace tqa = triqs::arrays;

using triqs::arrays::prod;
using triqs::arrays::sum;

#define MSG_PREFIX            __FILE__ << ":" << __LINE__ << ": "
#ifdef FK_MC_DEBUG
#define MY_DEBUG(MSG)            std::cout << std::boolalpha << MSG_PREFIX << MSG << std::endl; 
#else 
#define MY_DEBUG(MSG)
#endif
#define INFO(MSG)             std::cout << std::boolalpha << MSG << std::endl;
#define INFO2(MSG)            std::cout << "    " << std::boolalpha << MSG << std::endl;
#define INFO_NONEWLINE(MSG)   std::cout << MSG << std::flush;
#define ERROR(MSG)            std::cerr << MSG_PREFIX << MSG << std::endl;

typedef triqs::arrays::array<int,1> int_array_t;
typedef triqs::arrays::array<double,1> real_array_t;
typedef triqs::arrays::array_view<double,1> real_array_view_t;
typedef triqs::arrays::array<double,2> real_array2d_t;
typedef triqs::arrays::array_view<double,2> real_array2d_view_t;
typedef triqs::arrays::matrix<double>  real_matrix_t;
typedef triqs::arrays::matrix_view<double>  real_matrix_view_t;

template <typename T>
using EMatrixType = Eigen::Matrix<T,Eigen::Dynamic, Eigen::Dynamic, Eigen::AutoAlign|Eigen::RowMajor>;

typedef triqs::utility::parameters parameters;

template <size_t D>
inline std::ostream& operator<< (std::ostream& in, const std::array<size_t, D> arr){ 
    in << "{";  std::copy(arr.begin(), arr.end()-1, std::ostream_iterator<size_t>(in, ",")); in << *(arr.end()-1) << "}";
    return in;
};

/*
template <class arr_t>
inline double __prod(const arr_t& in){double out=1.0; out = std::accumulate(in.begin(), in.end(), 1.0, std::multiplies<double>()); return out; };
template <class arr_t>
inline double __sum(const arr_t& in){double out=0.0; out = std::accumulate(in.begin(), in.end(), 0.0, std::plus<double>()); return out; };
*/
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

namespace to_arr { 

template<int... Indices>
struct indices {
    using next = indices<Indices..., sizeof...(Indices)>;
};

template<int Size>
struct build_indices {
    using type = typename build_indices<Size - 1>::type::next;
};

template<>
struct build_indices<0> {
    using type = indices<>;
};

template<typename T>
using Bare = typename std::remove_cv<typename std::remove_reference<T>::type>::type;

template<typename Tuple>
constexpr
typename build_indices<std::tuple_size<Bare<Tuple>>::value>::type
make_indices()
{ return {}; }

template<typename Tuple, int... Indices>
std::array<
  typename std::tuple_element<0, Bare<Tuple>>::type,
    std::tuple_size<Bare<Tuple>>::value
>
to_array(Tuple&& tuple, indices<Indices...>)
{
    using std::get;
    return {{ get<Indices>(std::forward<Tuple>(tuple))... }};
}
} // end of namespace to_arr

template<typename Tuple>
auto to_array(Tuple&& tuple)
-> decltype( to_arr::to_array(std::declval<Tuple>(), to_arr::make_indices<Tuple>()) )
{
    return to_arr::to_array(std::forward<Tuple>(tuple), to_arr::make_indices<Tuple>());
}

}; // end of namespace FK

#endif // endif :: #ifndef __FK_MC_COMMON_H


