#pragma once
/** @file Direct.hpp
 * @brief Dispatch methods for P2P stage
 *
 */

#include "KernelTraits.hpp"
#include "executor/P2P.hpp"

#include <type_traits>

#include <vector>
#include <cassert>

class Direct {

public:

  /** Asymmetric matvec
   */
  template <typename Kernel,
            typename SourceIter, typename ChargeIter,
            typename TargetIter, typename ResultIter>
  inline static void matvec(const Kernel& K,
                            SourceIter s_first, SourceIter s_last,
                            ChargeIter c_first,
                            TargetIter t_first, TargetIter t_last,
                              ResultIter r_first)
  {
    P2P::block_eval(K.kernel(),
                    s_first, s_last, c_first,
                    t_first, t_last, r_first);
  }

  /** Symmetric matvec, off-diagonal block
   */
  template <typename Kernel,
            typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static void matvec(const Kernel& K,
                            SourceIter p1_first, SourceIter p1_last,
                            ChargeIter c1_first, ResultIter r1_first,
                            SourceIter p2_first, SourceIter p2_last,
                            ChargeIter c2_first, ResultIter r2_first)
  {
    P2P::block_eval(K.kernel(),
                    p1_first, p1_last, c1_first, r1_first,
                    p2_first, p2_last, c2_first, r2_first);
  }

  /** Symmetric matvec, diagonal block
   */
  template <typename Kernel,
            typename SourceIter, typename ChargeIter, typename ResultIter>
  inline static
  typename std::enable_if<!std::is_same<SourceIter,
                                        std::vector<typename Kernel::source_type>
                                        >::value
                          >::type   // XXX: Hack to avoid ambiguous declaration
  matvec(const Kernel& K,
         SourceIter p_first, SourceIter p_last,
         ChargeIter c_first, ResultIter r_first)
  {
    P2P::block_eval(K.kernel(),
                    p_first, p_last, c_first, r_first);
  }

  /** Convenience function for std::vector
   */
  template <typename Kernel>
  inline static void matvec(const Kernel& K,
                            const std::vector<typename Kernel::source_type>& s,
                            const std::vector<typename Kernel::charge_type>& c,
                            const std::vector<typename Kernel::target_type>& t,
                            std::vector<typename Kernel::result_type>& r)
  {
    assert(s.size() == c.size());
    assert(t.size() == r.size());

    // Pass to asymmetric p2p
    P2P::block_eval(K.kernel(),
                    s.begin(), s.end(), c.begin(),
                    t.begin(), t.end(), r.begin());
  }

  /** Convenience function for std::vector
   */
  template <typename Kernel>
  inline static void matvec(const Kernel& K,
                            const std::vector<typename Kernel::source_type>& p,
                            const std::vector<typename Kernel::charge_type>& c,
                            std::vector<typename Kernel::result_type>& r)
  {
    assert(p.size() == c.size());
    assert(p.size() == r.size());

    // Pass to symmetric p2p
    P2P::block_eval(K.kernel(),
                    p.begin(), p.end(), c.begin(), r.begin());
  }
};
