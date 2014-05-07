/** A useful suite for testing Kernel Expansions */

#include "fmmtl/dispatch/Dispatchers.hpp"
#include "fmmtl/Direct.hpp"

#include "LaplaceSpherical.hpp"
#include "YukawaCartesian.hpp"
//#include "StokesSpherical.hpp"

#include <iostream>


template <typename Expansion>
void two_level_test(const Expansion& K) {
  typedef Expansion expansion_type;
  typedef typename expansion_type::source_type source_type;
  typedef typename expansion_type::target_type target_type;
  typedef typename expansion_type::charge_type charge_type;
  typedef typename expansion_type::result_type result_type;
  typedef typename expansion_type::point_type point_type;
  typedef typename expansion_type::multipole_type multipole_type;
  typedef typename expansion_type::local_type local_type;

  // init source
  std::vector<source_type> s(1);
  s[0] = source_type(0, 0, 0);

  // init charge
  std::vector<charge_type> c(1);
  c[0] = charge_type(1);

  // init target
  std::vector<target_type> t(1);
  t[0] = target_type(0.98, 0.98, 0.98);

  // init results vectors for exact, FMM
  std::vector<result_type> rexact(1);
  rexact[0] = result_type(0);
  result_type rm2p = result_type(0);
  result_type rfmm = result_type(0);

  // test direct
  Direct::matvec(K, s, c, t, rexact);

  // setup initial multipole expansion
  multipole_type M;
  point_type M_center(0.05, 0.05, 0.05);
  point_type M_extent(0.1, 0.1, 0.1);
  INITM::apply(K, M, M_extent, 2u);
  K.S2M(s[0], c[0], M_center, M);

  // perform M2M
  multipole_type M2;
  point_type M2_center(0.1, 0.1, 0.1);
  point_type M2_extent(0.2, 0.2, 0.2);
  INITM::apply(K, M2, M2_extent, 1u);
  K.M2M(M, M2, M2_center - M_center);
  //K.S2M(s,c,M2_center,M2);

  // test M2L
  local_type L2;
  point_type L2_center(0.9, 0.9, 0.9);
  point_type L2_extent(0.2, 0.2, 0.2);
  INITL::apply(K, L2, L2_center, 1u);
  K.M2L(M2, L2, L2_center - M2_center);

  // test L2L
  local_type L;
  point_type L_center(0.95, 0.95, 0.95);
  point_type L_extent(0.1, 0.1, 0.1);
  INITL::apply(K, L, L_extent, 2u);
  K.L2L(L2, L, L_center - L2_center);

  // test L2T
  K.L2T(L2, L2_center, t[0], rfmm);

  // check errors
  std::cout << "rexact = " << rexact[0] << std::endl;
  std::cout << "rm2p = " << rm2p << "\n    "
            << "[" << (rm2p - rexact[0]) << "]" << std::endl;
  std::cout << "rfmm = " << rfmm << "\n    "
            << "[" << (rfmm - rexact[0]) << "]" << std::endl;
}


int main() {
  LaplaceSpherical K(5);
  //YukawaCartesian K(10, 0.1);
  //StokesSpherical K(5);

  two_level_test(K);
}
