

auto S2M = [&](int L, const source_box_type& sbox) {
  const point_type& s_center = sbox.center();
  const point_type s_scale = 1. / sbox.extents();

  // Precompute Lagrange interp matrix with sources transformed to ref grid
  auto make_ref = [&](const source_type& s) {
    return (s - s_center) * s_scale;
  };
  auto LgM = LagrangeMatrix<D,Q>(
      boost::make_transform_iterator(p_sources[sbox.body_begin()], make_ref),
      boost::make_transform_iterator(p_sources[sbox.body_end()],   make_ref));

  std::vector<complex_type> rhs;

  // For all the boxes in this level of the target tree
  auto tboxi = target_tree.box_begin(L);
  for (auto& M_AB : multipole[sbox]) {
    // The target center for this multipole
    const point_type& t_center = (*tboxi).center();
    ++tboxi;

    // Precompute the phase * charge
    rhs.assign(p_charges[sbox.body_begin()], p_charges[sbox.body_end()]);
    auto ri = std::begin(rhs);
    for (auto&& s : p_sources[sbox]) {
      *ri *= unit_polar(_M_ * K.phase(t_center, s));
      ++ri;
    }

    // Multiply  M_AB_i += LgM_ij rhs_j
    prod(LgM, rhs, M_AB);
  }
};


auto M2M = [&](int L, const source_box_type& sbox) {
  const point_type& s_center = sbox.center();
  const point_type  s_scale  = 1. / sbox.extents();

  for (source_box_type cbox : children(sbox)) {

    const point_type& c_center = cbox.center();
    const point_type& c_extent = cbox.extents();

    // Define the child box chebyshev grid
    std::array<point_type,pow_(Q,D)> c_cheb;
    auto cbi = c_cheb.begin();
    for (auto&& i : TensorIndexGridRange<D,Q>()) {
      for (unsigned d = 0; d < D; ++d)
        (*cbi)[d] = c_center[d] + Chebyshev<double,Q>::x[i[d]] * c_extent[d];
      ++cbi;
    }

    // Precompute Lagrange interp matrix with points transformed to ref grid
    // XXX: Avoid computing at all for quad/octrees
    auto make_ref = [&](const point_type& c) {
      return (c - s_center) * s_scale;
    };
    auto LgM = LagrangeMatrix<D,Q>(
        boost::make_transform_iterator(c_cheb.begin(), make_ref),
        boost::make_transform_iterator(c_cheb.end(),   make_ref));

    std::array<complex_type,pow_(Q,D)> rhs;

    // Accumulate into M_AB_t
    auto tboxi = target_tree.box_begin(L);
    for (auto&& M_AB : multipole[sbox]) {

      target_box_type tbox = *tboxi;
      ++tboxi;
      const point_type& t_center  = tbox.center();

      target_box_type pbox = tbox.parent();
      const point_type& p_center = pbox.center();
      int p_idx = pbox.index() - target_tree.box_begin(pbox.level())->index();

      // Precompute phase * M_ApBc   TODO: can reuse M_AB?
      rhs = multipole[cbox][p_idx];
      auto ri = std::begin(rhs);
      for (auto&& yi : c_cheb) {
        *ri *= unit_polar(_M_ * (K.phase(t_center, yi) - K.phase(p_center, yi)));
        ++ri;
      }

      // Multiply  M_AB_i += LgM_ij rhs_j
      prod(LgM, rhs, M_AB);
    }
  }
};


auto M2L = [&](int L) {
  // M2L and check the result
  int s_idx = -1;
  for (source_box_type sbox : boxes(L_max - L, source_tree)) {
    ++s_idx;

    const point_type& s_center = sbox.center();
    const point_type& s_extent = sbox.extents();

    // Define the source box Chebyshev grid
    std::array<point_type,pow_(Q,D)> s_cheb;
    {
    auto si = s_cheb.begin();
    for (auto&& i : TensorIndexGridRange<D,Q>()) {
      for (unsigned d = 0; d < D; ++d)
        (*si)[d] = s_center[d]
            + Chebyshev<double,Q>::x[i[d]] * s_extent[d];
      ++si;
    }
    }

    int t_idx = -1;
    for (target_box_type tbox : boxes(L, target_tree)) {
      ++t_idx;

      const point_type& t_center = tbox.center();
      const point_type& t_extent = tbox.extents();

      // Define the target box Chebyshev grid
      std::array<point_type,pow_(Q,D)> t_cheb;
      {
      auto ti = t_cheb.begin();
      for (auto&& i : TensorIndexGridRange<D,Q>()) {
        for (unsigned d = 0; d < D; ++d)
          (*ti)[d] = t_center[d]
              + Chebyshev<double,Q>::x[i[d]] * t_extent[d];
        ++ti;
      }
      }

      auto ti = t_cheb.begin();
      for (auto&& L_AB_i : local[tbox][s_idx]) {
        complex_type temp = 0;

        // Grab the multipole
        auto si = s_cheb.begin();
        for (auto&& M_AB_i : multipole[sbox][t_idx]) {
          temp += K.ampl(*ti,*si) * M_AB_i
              * unit_polar(_M_ * (K.phase(*ti,*si) -
                                  K.phase(t_center,*si)));
          ++si;
        }
        L_AB_i += temp * unit_polar(-_M_ * K.phase(*ti, s_center));
        ++ti;
      }
    }
  }
};


auto L2L = [&](int L, const target_box_type& tbox) {
  const point_type& t_center = tbox.center();
  const point_type  t_scale  = 1. / tbox.extents();

  for (target_box_type cbox : children(tbox)) {

    const point_type& c_center = cbox.center();
    const point_type& c_extent = cbox.extents();

    // Define the child box chebyshev grid
    std::array<point_type,pow_(Q,D)> c_cheb;
    auto cbi = c_cheb.begin();
    for (auto&& i : TensorIndexGridRange<D,Q>()) {
      for (unsigned d = 0; d < D; ++d)
        (*cbi)[d] = c_center[d] + Chebyshev<double,Q>::x[i[d]] * c_extent[d];
      ++cbi;
    }

    // Precompute Lagrange interp matrix with points transformed to ref grid
    // XXX: Avoid computing at all for quad/octrees
    auto make_ref = [&](const point_type& c) {
      return (c - t_center) * t_scale;
    };
    auto LgM = LagrangeMatrix<D,Q>(
        boost::make_transform_iterator(c_cheb.begin(), make_ref),
        boost::make_transform_iterator(c_cheb.end(),   make_ref));

    std::array<complex_type,pow_(Q,D)> lhs;

    // Accumulate from L_AB_t
    auto sboxi = source_tree.box_begin(L_max - L);
    for (auto&& L_AB : local[tbox]) {

      source_box_type sbox = *sboxi;
      ++sboxi;
      const point_type& s_center  = sbox.center();

      target_box_type pbox = sbox.parent();
      const point_type& p_center = pbox.center();
      int p_idx = pbox.index() - source_tree.box_begin(pbox.level())->index();

      // Multiply  lhs_j += LgM_ij L_AB_i
      lhs.fill(0);
      prod(trans(LgM), L_AB, lhs);

      // Precompute phase * L_AB   TODO: can reuse L_AB?
      auto li = std::begin(lhs);
      auto ci = std::begin(c_cheb);
      for (auto&& L_AcBp_i : local[cbox][p_idx]) {
        L_AcBp_i += *li *  unit_polar(_M_ * (K.phase(*ci, s_center) -
                                             K.phase(*ci, p_center)));
        ++li; ++ci;
      }
    }
  }
};



auto L2T = [&](int L, const target_box_type& tbox) {
  const point_type& t_center = tbox.center();
  const point_type  t_scale  = 1. / tbox.extents();

  // Precompute Lagrange interp matrix with targets transformed to ref grid
  auto make_ref = [&](const target_type& t) {
    return (t - t_center) * t_scale;
  };
  auto LgM = LagrangeMatrix<D,Q>(
      boost::make_transform_iterator(p_targets[tbox.body_begin()], make_ref),
      boost::make_transform_iterator(p_targets[tbox.body_end()],   make_ref));

  std::vector<result_type> temp(tbox.num_bodies());

  // For all the boxes in this level of the source tree
  auto sbi = source_tree.box_begin(L_max - L);
  for (const auto& L_AB : local[tbox]) {
    // The source center of L_AB
    const point_type& s_center = (*sbi).center();
    ++sbi;

    temp.assign(temp.size(), 0);

    // Multiply  temp_i += LgM_ij rhs_j
    prod(trans(LgM), L_AB, temp);

    auto ti = p_targets[tbox.body_begin()];
    auto tempi = std::begin(temp);
    for (auto&& ri : p_results[tbox]) {
      ri += *tempi * unit_polar(_M_ * (K.phase(*ti, s_center)));
      ++tempi;
      ++ti;
    }
  }
};
