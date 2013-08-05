#ifndef __FK_MC_COMMON_H
#define __FK_MC_COMMON_H

#include <triqs/arrays.hpp>
#include <triqs/arrays/algorithms.hpp>
#include <triqs/parameters/parameters.hpp>

#include <tuple>
#include <array>
#include <boost/math/special_functions/pow.hpp>
#include <numeric>

namespace fk {

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

}; // end of namespace FK

#endif // endif :: #ifndef __FK_MC_COMMON_H


