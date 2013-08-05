#ifndef __TRIQS_EXTRA_HPP
#define __TRIQS_EXTRA_HPP

#include <triqs/arrays.hpp>
#include <iostream>

/** All currently unsupported lib stuff. */
namespace triqs { namespace arrays { namespace extra { 

/** Weak storage. No ownership and no ref counting. */
template<typename ValueType>
struct weak_block {

    static constexpr bool is_weak = true;
    typedef ValueType value_type;
    ValueType* data_;
    size_t s_;
    explicit weak_block() { data_=nullptr; s_ = 0; }
    explicit weak_block(ValueType* in, size_t size):data_(in),s_(size) {}
    value_type & operator[](size_t i) const { return data_[i];}
    bool empty() const {return (s_ == 0);}
    size_t size() const {return s_;}

private:
    friend class boost::serialization::access;
    template<class Archive> void serialize(Archive & ar, const unsigned int version) { 
    //ar & boost::serialization::make_nvp("ptr",sptr); data_ = (sptr ? sptr->p : nullptr); s = (sptr ? sptr->size() : 0);
    for (size_t i=0; i<s_; i++) ar & data_[i]; // refer to boost
   }

};

 template <typename ValueType, int Rank, ull_t Opt=0, ull_t TraversalOrder= 0, bool Borrowed=false > class weak_view;

 // ---------------------- weak_view  --------------------------------

#define IMPL_TYPE indexmap_storage_pair< indexmaps::cuboid::map<Rank,Opt,TraversalOrder>, \
 weak_block<ValueType>, Opt, TraversalOrder, Tag::array_view > 

 template <typename ValueType, int Rank, ull_t Opt, ull_t TraversalOrder, bool Borrowed>
  class weak_view : Tag::array_view, TRIQS_MODEL_CONCEPT(MutableCuboidArray), public IMPL_TYPE {
   static_assert( Rank>0, " Rank must be >0");
   public:   
   typedef typename IMPL_TYPE::indexmap_type indexmap_type;
   typedef typename IMPL_TYPE::storage_type storage_type;
   typedef array     <ValueType,Rank,Opt,TraversalOrder>       non_view_type;
   typedef weak_view<ValueType,Rank,Opt,TraversalOrder>       view_type;
   typedef weak_view<ValueType,Rank,Opt,TraversalOrder,true>  weak_view_type;
   typedef void has_view_type_tag;
 
   /// Build from an IndexMap and a storage 
   template<typename S> weak_view (indexmap_type const & Ind,S const & Mem): IMPL_TYPE(Ind, Mem) {}

   /// Copy constructor
   weak_view(weak_view const & X): IMPL_TYPE(X.indexmap(),X.storage()) {}

   /// Build from anything that has an indexmap and a storage compatible with this class
   template<typename ISP> weak_view(const ISP & X): IMPL_TYPE(X.indexmap(),X.storage()) {}

#ifdef TRIQS_WITH_PYTHON_SUPPORT
   /// Build from a numpy.array (X is a borrowed reference) : throws if X is not a numpy.array 
   explicit weak_view (PyObject * X): IMPL_TYPE(X, false, "weak_view "){}
#endif

   weak_view () = delete;

   // Move
   weak_view(weak_view && X) { this->swap_me(X); }

   /// Swap
   friend void swap( weak_view & A, weak_view & B) { A.swap_me(B);}

   /// Rebind the view
   void rebind (weak_view const & X) { this->indexmap_ = X.indexmap_; this->storage_ = X.storage_;}

   /// Assignment. The size of the array MUST match exactly, except in the empty case 
   template<typename RHS> weak_view & operator=(RHS const & X) { triqs_arrays_assign_delegation(*this,X); return *this; }

   ///
   weak_view & operator=(weak_view const & X) { triqs_arrays_assign_delegation(*this,X); return *this; } //without this, the standard = is synthetized...

   // Move assignment not defined : will use the copy = since view must copy data

   TRIQS_DEFINE_COMPOUND_OPERATORS(weak_view);
   // to forbid serialization of views...
   //template<class Archive> void serialize(Archive & ar, const unsigned int version) = delete;
  
  };

#undef IMPL_TYPE

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
 
}}}

namespace fk {
    typedef triqs::arrays::extra::weak_view<double, 1> real_array_weak_view_t;
    typedef triqs::arrays::extra::weak_view<double, 2> real_array2d_weak_view_t;
}

#endif // endif :: #ifndef __TRIQS_EXTRA_HPP
