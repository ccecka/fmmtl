#pragma once

#include <boost/iterator/iterator_adaptor.hpp>

#include <iostream>

#include "fmmtl/meta/func_traits.hpp"


template <typename Kernel>
struct KernelTraits {
 public:
  typedef KernelTraits<Kernel>                    self_type;
  typedef typename Kernel::kernel_type            kernel_type;

  typedef typename kernel_type::source_type       source_type;
  typedef typename kernel_type::target_type       target_type;
  typedef typename kernel_type::charge_type       charge_type;
  typedef typename kernel_type::kernel_value_type kernel_value_type;
  typedef typename kernel_type::result_type       result_type;

  // Kernel evaluation operator, K(t,s)
  HAS_MEM_FUNC(HasEvalOp,
               kernel_value_type, operator(),
               const target_type&, const source_type&);
  static const bool has_eval_op = HasEvalOp<Kernel>::value;
  // Kernel transpose, K.transpose(K(t,s))
  HAS_MEM_FUNC(HasTranspose,
               kernel_value_type, transpose,
               const kernel_value_type&);
  static const bool has_transpose = HasTranspose<Kernel>::value;

 protected:
    // A dummy iterator adaptor to check for templated vectorized methods
  template <typename T>
  struct dumb_iterator : public boost::iterator_adaptor<dumb_iterator<T>,T*> {};

  typedef dumb_iterator<source_type> source_iterator;
  typedef dumb_iterator<charge_type> charge_iterator;
  typedef dumb_iterator<target_type> target_iterator;
  typedef dumb_iterator<result_type> result_iterator;

 public:
  HAS_MEM_FUNC(HasP2Psymm,
               void, P2P,
               source_iterator, source_iterator, charge_iterator,
               target_iterator, target_iterator, charge_iterator,
               result_iterator, result_iterator);
  static const bool has_vector_P2P_symm = HasP2Psymm<Kernel>::value;
  HAS_MEM_FUNC(HasP2Pasymm,
               void, P2P,
               source_iterator, source_iterator, charge_iterator,
               target_iterator, target_iterator, result_iterator);
  static const bool has_vector_P2P_asymm = HasP2Pasymm<Kernel>::value;

  friend std::ostream& operator<<(std::ostream& s, const self_type& traits) {
    s << "has_eval_op: "           << traits.has_eval_op           << std::endl;
    s << "has_transpose: "         << traits.has_transpose         << std::endl;
    s << "has_vector_P2P_symm: "   << traits.has_vector_P2P_symm   << std::endl;
    s << "has_vector_P2P_asymm: "  << traits.has_vector_P2P_asymm;
    return s;
  }
};


#include "fmmtl/meta/dimension.hpp"

template <typename Expansion>
struct ExpansionTraits
    : public KernelTraits<typename Expansion::kernel_type> {
  typedef Expansion                               expansion_type;
  typedef typename expansion_type::kernel_type    kernel_type;
  typedef ExpansionTraits<expansion_type>         self_type;
  typedef KernelTraits<kernel_type>               super_type;

  typedef typename super_type::source_type        source_type;
  typedef typename super_type::target_type        target_type;
  typedef typename super_type::charge_type        charge_type;
  typedef typename super_type::kernel_value_type  kernel_value_type;
  typedef typename super_type::result_type        result_type;

  typedef typename super_type::source_iterator    source_iterator;
  typedef typename super_type::charge_iterator    charge_iterator;
  typedef typename super_type::target_iterator    target_iterator;
  typedef typename super_type::result_iterator    result_iterator;

  typedef typename expansion_type::point_type     point_type;
  static const std::size_t dimension = fmmtl::dimension<point_type>::value;

  typedef typename expansion_type::multipole_type multipole_type;
  typedef typename expansion_type::local_type     local_type;

  // Converters
  HAS_MEM_FUNC(HasSourcePoint,
               point_type, source_point,
               const source_type&);
  static const bool has_source_point = HasSourcePoint<Expansion>::value;
  HAS_MEM_FUNC(HasTargetPoint,
               point_type, target_point,
               const target_type&);
  static const bool has_target_point = HasTargetPoint<Expansion>::value;

  // Initializers
  HAS_MEM_FUNC(HasInitMultipole,
               void, init_multipole,
               multipole_type&, const point_type&, unsigned);
  static const bool has_init_multipole = HasInitMultipole<Expansion>::value;
  HAS_MEM_FUNC(HasInitLocal,
               void, init_local,
               local_type&, const point_type&, unsigned);
  static const bool has_init_local = HasInitLocal<Expansion>::value;

  // P2M
  HAS_MEM_FUNC(HasScalarP2M,
               void, P2M,
               const source_type&, const charge_type&,
               const point_type&, multipole_type&);
  static const bool has_scalar_P2M = HasScalarP2M<Expansion>::value;
  HAS_MEM_FUNC(HasVectorP2M,
               void, P2M,
               source_iterator, source_iterator, charge_iterator,
               const point_type&, multipole_type&);
  static const bool has_vector_P2M = HasVectorP2M<Expansion>::value;
  static const bool has_P2M = (has_scalar_P2M || has_vector_P2M);

  // P2L
  HAS_MEM_FUNC(HasScalarP2L,
               void, P2L,
               const source_type&, const charge_type&,
               const point_type&, local_type&);
  static const bool has_scalar_P2L = HasScalarP2L<Expansion>::value;
  HAS_MEM_FUNC(HasVectorP2L,
               void, P2L,
               source_iterator, source_iterator, charge_iterator,
               const point_type&, local_type&);
  static const bool has_vector_P2L = HasVectorP2L<Expansion>::value;
  static const bool has_P2L = (has_scalar_P2L || has_vector_P2L);

  // M2M
  HAS_MEM_FUNC(HasM2M,
               void, M2M,
               const multipole_type&, multipole_type&, const point_type&);
  static const bool has_M2M = HasM2M<Expansion>::value;

  // M2L
  HAS_MEM_FUNC(HasM2L,
               void, M2L,
               const multipole_type&, local_type&, const point_type&);
  static const bool has_M2L = HasM2L<Expansion>::value;

  // MAC
  HAS_MEM_FUNC(HasDynMAC,
               bool, MAC,
               const multipole_type&, const local_type&);
  static const bool has_dynamic_MAC = HasDynMAC<Expansion>::value;

  // L2L
  HAS_MEM_FUNC(HasL2L,
               void, L2L,
               const local_type&, local_type&, const point_type&);
  static const bool has_L2L = HasL2L<Expansion>::value;

  // M2P
  HAS_MEM_FUNC(HasScalarM2P,
               void, M2P,
               const multipole_type&, const point_type&,
               const target_type&, result_type&);
  static const bool has_scalar_M2P = HasScalarM2P<Expansion>::value;
  HAS_MEM_FUNC(HasVectorM2P,
               void, M2P,
               const multipole_type&, const point_type&,
               target_iterator, target_iterator, result_iterator);
  static const bool has_vector_M2P = HasVectorM2P<Expansion>::value;
  static const bool has_M2P = (has_scalar_M2P || has_vector_M2P);

  // L2P
  HAS_MEM_FUNC(HasScalarL2P,
               void, L2P,
               const local_type&, const point_type&,
               const target_type&, result_type&);
  static const bool has_scalar_L2P = HasScalarL2P<Expansion>::value;
  HAS_MEM_FUNC(HasVectorL2P,
               void, L2P,
               const local_type&, const point_type&,
               target_iterator, target_iterator, result_iterator);
  static const bool has_vector_L2P = HasVectorL2P<Expansion>::value;
  static const bool has_L2P = (has_scalar_L2P || has_vector_L2P);

  friend std::ostream& operator<<(std::ostream& s, const self_type& traits) {
    s << static_cast<super_type>(traits)                     << std::endl;
    s << "has_init_multipole: " << traits.has_init_multipole << std::endl;
    s << "has_init_local: "     << traits.has_init_local     << std::endl;
    s << "has_P2M: "            << traits.has_P2M            << std::endl;
    s << "  has_scalar_P2M: "   << traits.has_scalar_P2M     << std::endl;
    s << "  has_vector_P2M: "   << traits.has_vector_P2M     << std::endl;
    s << "has_P2L: "            << traits.has_P2L            << std::endl;
    s << "  has_scalar_P2L: "   << traits.has_scalar_P2L     << std::endl;
    s << "  has_vector_P2L: "   << traits.has_vector_P2L     << std::endl;
    s << "has_M2M: "            << traits.has_M2M            << std::endl;
    s << "has_M2L: "            << traits.has_M2L            << std::endl;
    s << "has_L2L: "            << traits.has_L2L            << std::endl;
    s << "has_M2P: "            << traits.has_M2P            << std::endl;
    s << "  has_scalar_M2P: "   << traits.has_scalar_M2P     << std::endl;
    s << "  has_vector_M2P: "   << traits.has_vector_M2P     << std::endl;
    s << "has_L2P: "            << traits.has_M2P            << std::endl;
    s << "  has_scalar_L2P: "   << traits.has_scalar_L2P     << std::endl;
    s << "  has_vector_L2P: "   << traits.has_vector_L2P     << std::endl;
    s << "has_dynamic_MAC: "    << traits.has_dynamic_MAC;
    return s;
  }
};

#define FMMTL_IMPORT_KERNEL_TRAITS(K)                                     \
  typedef typename KernelTraits<K>::kernel_type        kernel_type;       \
  typedef typename KernelTraits<K>::kernel_value_type  kernel_value_type; \
  typedef typename KernelTraits<K>::source_type        source_type;       \
  typedef typename KernelTraits<K>::target_type        target_type;       \
  typedef typename KernelTraits<K>::charge_type        charge_type;       \
  typedef typename KernelTraits<K>::result_type        result_type

#define FMMTL_IMPORT_EXPANSION_TRAITS(E)                                  \
  FMMTL_IMPORT_KERNEL_TRAITS(E);                                          \
  typedef typename ExpansionTraits<E>::expansion_type     expansion_type; \
  typedef typename ExpansionTraits<E>::multipole_type     multipole_type; \
  typedef typename ExpansionTraits<E>::local_type         local_type;     \
  typedef typename ExpansionTraits<E>::point_type         point_type
