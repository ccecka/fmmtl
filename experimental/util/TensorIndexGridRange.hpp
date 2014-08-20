#pragma once


/** Quickie class for iterating over the integer tensor range
 * (0,0,...,0) x (Q,Q,...,Q)    (where cardinality is DIM)
 */
template <std::size_t DIM, std::size_t Q>
struct TensorIndexGridRange {
  typedef std::array<unsigned,DIM> value_type;
  value_type i_ = {{}};

  // Prevent copying the value_type when iterating
  struct Reference {
    TensorIndexGridRange<DIM, Q>& t;

    // Increment the grid index and carry into the next dimensions
    void operator++() {
      for (std::size_t k = 0; ++t.i_[k] == Q && k != DIM-1; ++k)
        t.i_[k] = 0;
    }
    // Current != end of the range
    template <typename T>
    bool operator!=(T&&) const {
      return t.i_[DIM-1] != Q;
    }
    // Return the grid index
    const value_type& operator*() const {
      return t.i_;
    }
  };

  Reference begin() { return Reference{*this}; }
  Reference end()   { return Reference{*this}; }
};
