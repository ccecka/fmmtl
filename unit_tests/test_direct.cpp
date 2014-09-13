#include "fmmtl/Kernel.hpp"
#include "fmmtl/Direct.hpp"

#include <vector>
#include <iostream>

struct TempKernel
    : public fmmtl::Kernel<TempKernel> {
  //! Source type
  typedef double source_type;
  //! Charge type
  typedef double charge_type;
  //! Target type
  typedef double target_type;
  //! Result type
  typedef double result_type;

  //! The return type of a kernel evaluation
  typedef unsigned kernel_value_type;

  kernel_value_type operator()(const target_type& t,
                               const source_type& s) const {
    (void) t;
    (void) s;

    std::cout << "In op()\n";

    return kernel_value_type(1);
  }

  kernel_value_type transpose(const kernel_value_type& kst) const {
    (void) kst;
    std::cout << "In Transpose\n";
    return kernel_value_type(1);
  }

  // A fake S2T that shouldn't be called
  void S2T(int a) const {
    (void) a;
    std::cout << "In S2T(int)\n";
  }

#if 0
  template <typename PointIter, typename ChargeIter, typename ResultIter>
  void S2T(PointIter s_begin, PointIter s_end, ChargeIter c_begin,
           PointIter t_begin, PointIter t_end, ResultIter r_begin) const {
    (void) s_begin;
    (void) s_end;
    (void) c_begin;
    (void) t_begin;
    (void) t_end;
    (void) r_begin;

    std::cout << "In vector S2T\n";
  }
#endif
};




int main() {
  typedef TempKernel kernel_type;
  typedef kernel_type::source_type source_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::target_type target_type;
  typedef kernel_type::result_type result_type;

  kernel_type K;

  // init source
  std::vector<source_type> sources(2);

  // init charge
  std::vector<charge_type> charges(sources.size());

  // init target
  std::vector<target_type> targets(2);

  // init results vectors for exact
  std::vector<result_type> exact(targets.size());

  // test direct
  fmmtl::direct(K, sources, charges, targets, exact);
  std::cout << exact[0] << "\n";
  fmmtl::direct(K, sources, charges, exact);
  std::cout << exact[0] << "\n";
}
