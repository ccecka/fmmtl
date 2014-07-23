/** A useful suite for testing Kernel Expansions */

#include "fmmtl/dispatch/Dispatchers.hpp"
#include "fmmtl/Direct.hpp"

#include "LaplaceSpherical.hpp"
#include "YukawaCartesian.hpp"
//#include "StokesSpherical.hpp"

#include <iostream>

template <typename Expansion>
void single_level_test(const Expansion& K) {
  typedef Expansion expansion_type;
  typedef typename expansion_type::point_type point_type;
  typedef typename expansion_type::source_type source_type;
  typedef typename expansion_type::target_type target_type;
  typedef typename expansion_type::charge_type charge_type;
  typedef typename expansion_type::result_type result_type;
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
  point_type M_center(0.1, 0.1, 0.1);
  point_type M_extent(0.2, 0.2, 0.2);
  INITM::apply(K, M, M_extent, 1u);
  K.S2M(s[0], c[0], M_center, M);

  // test M2T
  //K.M2T(M, M_center, t[0], rm2p);

  // test M2L, L2T
  local_type L;
  point_type L_center(0.9, 0.9, 0.9);
  point_type L_extent(0.2, 0.2, 0.2);
  auto d = L_center - M_center;
  printf("DIST: (%lg, %lg, %lg) : %lg\n", d[0], d[1], d[2], norm_2(d));
  INITL::apply(K, L, L_extent, 1u);
  K.M2L(M, L, L_center - M_center);
  K.L2T(L, L_center, t[0], rfmm);

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

  single_level_test(K);
}
