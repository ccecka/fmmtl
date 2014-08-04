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

  // For all the boxes in this level of the target tree
  auto tboxi = target_tree.box_begin(L);
  for (auto& M_AB : multipole[sbox]) {
    // The target center for this multipole
    const point_type& t_center = (*tboxi).center();
    ++tboxi;

    // TODO: Precompute the phase * charge

    // For each multiindex
    auto mabi = M_AB.begin();
    for (auto&& i : TensorIndexGridRange<D,Q>()) {

      // Accumulate into M_AB_i
      complex_type& M_AB_i = *mabi;
      ++mabi;

      auto ci = p_charges[sbox.body_begin()];
      std::size_t j = 0;   // XXX, Abstract as mat-vec
      for (auto&& s : p_sources[sbox]) {
        M_AB_i += unit_polar(_M_ * K.phase(t_center, s)) * LgM(i,j) * (*ci);
        ++ci; ++j;
      }
    }
  }
};


auto M2M = [&](int L, const source_box_type& cbox) {
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

  source_box_type sbox = cbox.parent();

  const point_type& s_center = sbox.center();
  const point_type  s_scale  = 1. / sbox.extents();

  // Precompute Lagrange interp matrix with points transformed to ref grid
  // XXX: Avoid computing at all for quad/octrees
  auto make_ref = [&](const point_type& c) {
    return (c - s_center) * s_scale;
  };
  auto LgM = LagrangeMatrix<D,Q>(
      boost::make_transform_iterator(c_cheb.begin(), make_ref),
      boost::make_transform_iterator(c_cheb.end(),   make_ref));

  // Accumulate into M_AB_t
  auto tboxi = target_tree.box_begin(L);
  for (auto&& M_AB : multipole[sbox]) {

    target_box_type tbox = *tboxi;
    ++tboxi;
    const point_type& t_center  = tbox.center();

    target_box_type pbox = tbox.parent();
    const point_type& p_center = pbox.center();
    int p_idx = pbox.index() - target_tree.box_begin(pbox.level())->index();

    // Accumulate
    auto mabi = M_AB.begin();
    for (auto&& i : TensorIndexGridRange<D,Q>()) {

      // Accumulate into M_AB_i
      complex_type& M_AB_i = *mabi;
      ++mabi;

      // For each element of i_prime
      auto yi = c_cheb.begin();
      std::size_t j = 0;  // Lift the matvec
      for (auto&& M_ApBc_ip : multipole[cbox][p_idx]) {
        M_AB_i += M_ApBc_ip  * LgM(i,j)
            * unit_polar(_M_ * (K.phase(t_center, *yi) -
                                K.phase(p_center, *yi)));
        ++yi; ++j;
      }
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

        // Grab the multipole
        auto si = s_cheb.begin();
        for (auto&& M_AB_i : multipole[sbox][t_idx]) {
          L_AB_i += K.ampl(*ti,*si) * M_AB_i
              * unit_polar(_M_ * (K.phase(*ti,*si) - K.phase(t_center,*si)));
          ++si;
        }
        L_AB_i *= unit_polar(-_M_ * K.phase(*ti, s_center));
        ++ti;
      }
    }
  }
};


auto L2L = [&](int L, const target_box_type& tbox) {
  const point_type& t_center = tbox.center();
  const point_type& t_extent = tbox.extents();

  // Define the target box Chebyshev grid
  std::array<point_type,pow_(Q,D)> t_cheb;
  auto cbi = t_cheb.begin();
  for (auto&& i : TensorIndexGridRange<D,Q>()) {
    for (unsigned d = 0; d < D; ++d)
      (*cbi)[d] = t_center[d]
          + Chebyshev<double,Q>::x[i[d]] * t_extent[d];
    ++cbi;
  }

  // Get the target box's parent
  target_box_type pbox = tbox.parent();

  const point_type& p_center = pbox.center();
  const point_type  p_scale  = 1. / pbox.extents();

  // Precompute Lagrange interp matrix with points transformed to ref grid
  // XXX: Avoid computing at all for quad/octrees
  auto make_ref = [&](const point_type& c) {
    return (c - p_center) * p_scale;
  };
  auto LgM = LagrangeMatrix<D,Q>(
      boost::make_transform_iterator(t_cheb.begin(), make_ref),
      boost::make_transform_iterator(t_cheb.end(),   make_ref));

  // For all the child boxes
  int c_idx = -1;
  for (source_box_type cbox : boxes(L_max - L + 1, source_tree)) {
    ++c_idx;

    const point_type& c_center = cbox.center();

    // Get the parent
    source_box_type sbox = cbox.parent();
    int s_idx = sbox.index() - boxes(L_max - L, source_tree).begin()->index();
    const point_type& s_center = sbox.center();

    // Accumulate
    int i_idx = -1;
    for (auto&& i : TensorIndexGridRange<D,Q>()) {
      ++i_idx;

      complex_type& L_ApBc_ip = local[pbox][c_idx][i_idx];

      // For each element of i_prime
      auto xi = t_cheb.begin();
      std::size_t j = 0;  // Lift the matvec
      for (auto&& L_AB_i : local[tbox][s_idx]) {
        L_AB_i += L_ApBc_ip * LgM(i,j)
            * unit_polar(_M_ * (K.phase(*xi, c_center) -
                                K.phase(*xi, s_center)));
        ++xi;
        ++j;
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

  // For all the boxes in this level of the source tree
  auto sbi = source_tree.box_begin(L_max - L);
  for (const auto& L_AB : local[tbox]) {
    // The source center of L_AB
    const point_type& s_center = (*sbi).center();
    ++sbi;

    // For each multiindex
    auto labi = L_AB.begin();
    for (auto&& i : TensorIndexGridRange<D,Q>()) {

      // Compute from L_AB_i
      const complex_type& L_AB_i = *labi;
      ++labi;

      auto ri = p_results[tbox.body_begin()];
      std::size_t j = 0;     // XXX, Abstract as mat-vec
      for (auto&& t : p_targets[tbox]) {
        *ri += LgM(i,j) * L_AB_i * unit_polar(_M_ * K.phase(t, s_center));
        ++j; ++ri;
      }
    }
  }
};
