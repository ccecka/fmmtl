#pragma once

#include <boost/iterator/iterator_adaptor.hpp>

#include <iostream>
#include <type_traits>
#include <utility>

// TODO: Make C++03 compatable?
#define SFINAE_TEMPLATE(NAME, KERNEL, OP)                               \
  template <typename ReturnType, typename... Arg>                       \
  struct NAME {                                                         \
    template <class A, class S = void>                                  \
    struct has_op : std::false_type {};                                 \
    template <class A>                                                  \
    struct has_op<A,                                                    \
      typename std::enable_if<                                          \
        std::is_same<                                                   \
          ReturnType,                                                   \
          decltype(std::declval<const A&>().OP(std::declval<Arg>()...)) \
        >::value                                                        \
      >::type> : std::true_type {};                                     \
    static constexpr bool value = has_op<KERNEL>::value;                \
  }


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

  // Kernel Evaluations and P2P
  SFINAE_TEMPLATE(HasEvalOp,kernel_type,operator());
  static constexpr bool has_eval_op =
      HasEvalOp<kernel_value_type,
                const target_type&, const source_type&>::value;
  SFINAE_TEMPLATE(HasTranspose,kernel_type,transpose);
  static constexpr bool has_transpose =
      HasTranspose<kernel_value_type,
                   const kernel_value_type&>::value;

 protected:
    // A dummy iterator adaptor to check for templated vectorized methods
  template <typename T>
  struct dumb_iterator : public boost::iterator_adaptor<dumb_iterator<T>,T*> {};

  typedef dumb_iterator<source_type> source_iterator;
  typedef dumb_iterator<charge_type> charge_iterator;
  typedef dumb_iterator<target_type> target_iterator;
  typedef dumb_iterator<result_type> result_iterator;

 public:
  SFINAE_TEMPLATE(HasP2P,kernel_type,P2P);
  static constexpr bool has_vector_P2P_symm =
      HasP2P<void,
             source_iterator, source_iterator, charge_iterator,
             target_iterator, target_iterator, charge_iterator,
             result_iterator, result_iterator>::value;
  static constexpr bool has_vector_P2P_asymm =
      HasP2P<void,
             source_iterator, source_iterator, charge_iterator,
             target_iterator, target_iterator, result_iterator>::value;

  friend std::ostream& operator<<(std::ostream& s, const self_type& traits) {
    s << "has_eval_op: "           << traits.has_eval_op           << std::endl;
    s << "has_transpose: "         << traits.has_transpose         << std::endl;
    s << "has_vector_P2P_symm: "   << traits.has_vector_P2P_symm   << std::endl;
    s << "has_vector_P2P_asymm: "  << traits.has_vector_P2P_asymm;
    return s;
  }
};




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

  static constexpr unsigned dimension = expansion_type::dimension;
  typedef typename expansion_type::point_type     point_type;
  // TODO: Check point_type = Vec<dimension,double> or generalize

  typedef typename expansion_type::multipole_type multipole_type;
  typedef typename expansion_type::local_type     local_type;

  // Initializers
  SFINAE_TEMPLATE(HasInitMultipole,expansion_type,init_multipole);
  static constexpr bool has_init_multipole =
      HasInitMultipole<void,
                       multipole_type&, const point_type&, unsigned>::value;
  SFINAE_TEMPLATE(HasInitLocal,expansion_type,init_local);
  static constexpr bool has_init_local =
      HasInitLocal<void,
                   local_type&, const point_type&, unsigned>::value;

  // P2M
  SFINAE_TEMPLATE(HasP2M,expansion_type,P2M);
  static constexpr bool has_scalar_P2M =
      HasP2M<void,
             const source_type&, const charge_type&,
             const point_type&, multipole_type&>::value;
  static constexpr bool has_vector_P2M =
      HasP2M<void,
             source_iterator, source_iterator, charge_iterator,
             const point_type&, multipole_type&>::value;
  static constexpr bool has_P2M =
      (has_scalar_P2M || has_vector_P2M);

  // P2L
  SFINAE_TEMPLATE(HasP2L,expansion_type,P2M);
  static constexpr bool has_scalar_P2L =
      HasP2L<void,
             const source_type&, const charge_type&,
             const point_type&, local_type&>::value;
  static constexpr bool has_vector_P2L =
      HasP2L<void,
             source_iterator, source_iterator, charge_iterator,
             const point_type&, local_type&>::value;
  static constexpr bool has_P2L =
      (has_scalar_P2L || has_vector_P2L);

  // M2M
  SFINAE_TEMPLATE(HasM2M,expansion_type,M2M);
  static constexpr bool has_M2M =
      HasM2M<void,
             const multipole_type&, multipole_type&, const point_type&>::value;

  // M2L
  SFINAE_TEMPLATE(HasM2L,expansion_type,M2L);
  static constexpr bool has_M2L =
      HasM2L<void,
             const multipole_type&, local_type&, const point_type&>::value;

  // MAC
  SFINAE_TEMPLATE(HasDynMAC,expansion_type,MAC);
  static constexpr bool has_dynamic_MAC =
      HasDynMAC<bool,
                const multipole_type&, const local_type&>::value;

  // L2L
  SFINAE_TEMPLATE(HasL2L,expansion_type,L2L);
  static constexpr bool has_L2L =
      HasL2L<void,
             const local_type&, local_type&, const point_type&>::value;

  // M2P
  SFINAE_TEMPLATE(HasM2P,expansion_type,M2P);
  static constexpr bool has_scalar_M2P =
      HasM2P<void,
             const multipole_type&, const point_type&,
             const target_type&, result_type&>::value;
  static constexpr bool has_vector_M2P =
      HasM2P<void,
             const multipole_type&, const point_type&,
             target_iterator, target_iterator, result_iterator>::value;
  static constexpr bool has_M2P =
      (has_scalar_M2P || has_vector_M2P);

  // L2P
  SFINAE_TEMPLATE(HasL2P,expansion_type,L2P);
  static constexpr bool has_scalar_L2P =
      HasL2P<void,
             const local_type&, const point_type&,
             const target_type&, result_type&>::value;
  static constexpr bool has_vector_L2P =
      HasL2P<void,
             const local_type&, const point_type&,
             target_iterator, target_iterator, result_iterator>::value;
  static constexpr bool has_L2P =
      (has_scalar_L2P || has_vector_L2P);

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

#undef SFINAE_TEMPLATE

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

