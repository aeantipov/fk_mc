#ifndef __FK_MC_COMMON_H
#define __FK_MC_COMMON_H

#include <triqs/arrays.hpp>
#include <triqs/parameters/parameters.hpp>

#include <tuple>
#include <array>
#include <boost/math/special_functions/pow.hpp>

namespace fk {

#define MSG_PREFIX            __FILE__ << ":" << __LINE__ << ": "
#define DEBUG(MSG)            std::cout << std::boolalpha << MSG_PREFIX << MSG << std::endl;
#define INFO(MSG)             std::cout << std::boolalpha << MSG << std::endl;
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

}; // end of namespace FK

#endif // endif :: #ifndef __FK_MC_COMMON_H


