#ifndef __TRIQS_EXTRA_HPP
#define __TRIQS_EXTRA_HPP

#include <triqs/arrays.hpp>
#include <iostream>

/** All currently unsupported lib stuff. */
namespace triqs { namespace arrays { 

 // reinterpret_array from Olivier 
 template <typename V, int R, ull_t Opt, ull_t To, typename ... I>  
  array_view<V, sizeof...(I), Opt> reinterpret_array_view (array<V,R,Opt,To> const & a, I ... index) { 
   static int constexpr rank = sizeof...(I);
   typedef array_view<V, rank, Opt, indexmaps::mem_layout::c_order(rank)> return_t;
   return return_t (typename return_t::indexmap_type (mini_vector<size_t,rank>(index...)) , a.storage());
  }


 template<typename T> class immutable_diagonal_matrix_view : TRIQS_MODEL_CONCEPT(ImmutableMatrix) { 
  array_view<T,1> data;
  public:

  //immutable_diagonal_matrix_view(size_t i=1) : data ( array<T,1>(i)) {}
  immutable_diagonal_matrix_view(array_view<T,1> v) : data (v) {}
 
  // the ImmutableMatrix concept 
  typedef T value_type;
  T operator()(size_t i, size_t j) const { return (i==j ? data(i) : 0);}

  typedef indexmaps::cuboid::domain_t<2> domain_type;
  domain_type domain() const { auto s = data.shape()[0]; return mini_vector<size_t,2>(s,s);}

  //
  friend std::ostream & operator<<(std::ostream & out, immutable_diagonal_matrix_view const & d) { return out<<"diagonal_matrix "<<d.data;}
  
  // ----------------------
  // should be remove from concept. redundant....
  // need to clean dim0, dim1 and shape and make them free function everywhere (deduced from domain)
  size_t dim0() const { return data.shape()[0];}
  size_t dim1() const { return data.shape()[0];}
  mini_vector<size_t,2> shape() const { auto s = data.shape()[0]; return mini_vector<size_t,2>(s,s);}
 
 };
 
}}



#endif // endif :: #ifndef __TRIQS_EXTRA_HPP
