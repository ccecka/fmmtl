#pragma once

/* @class      : HAS_MEM_FUNC
 * @brief      : This macro may be used to check if a class has a public
 *               const member function with particular signature.
 * @param NAME        : Name of struct this macro creates.
 * @param RETURN_TYPE : Return type of the member function to test for.
 * @param FUNC        : Name of member function to test for.
 * @param (...)       : The argument types of the member function to test for.
 *                      These complete the signature of the member funtion.
 * @note This macro is C++03 compatible
 */
#define HAS_MEM_FUNC(NAME, RETURN_TYPE, FUNC, ...)                      \
  template <typename CLASS>                                             \
  struct NAME {                                                         \
    typedef RETURN_TYPE (CLASS::*A)(__VA_ARGS__) const;                 \
    template <typename U, U> struct type_check;                         \
    template <typename U> static char chk(type_check<A, &U::FUNC>*);    \
    template <typename  > static int  chk(...);                         \
    static bool const value = sizeof(chk<CLASS>(0)) == sizeof(char);    \
  }



/** @brief : Another version that uses C++11
#include <type_traits>
#include <utility>

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
*/
