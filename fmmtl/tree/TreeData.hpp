#pragma once


/** Maps boxes and box iterators to data and data iterators
 */
template <typename T, typename Tree>
struct BoxBind {
  typedef std::vector<T> container_type;
  typedef typename container_type::iterator iterator;
  typedef typename container_type::const_iterator const_iterator;

  container_type data;

  BoxBind(const Tree& tree, unsigned size = tree.boxes())
      : data(size) {
  }

  T& operator[](const typename Tree::box_type& box) {
    return data[box.index()];
  }
  const T& operator[](const typename Tree::box_type& box) const {
    return data[box.index()];
  }

  iterator operator[](const typename Tree::box_iterator& bi) {
    return data.begin() + (*bi).index();
  }
  const_iterator operator[](const typename Tree::box_iterator& bi) const {
    return data.begin() + (*bi).index();
  }
};




/** Maps bodies and body iterators to data and data iterators
 */
template <typename T, typename Tree>
struct BodyBind {
  typedef std::vector<T> container_type;
  typedef typename container_type::iterator iterator;
  typedef typename container_type::const_iterator const_iterator;

  const Tree& tree_;
  container_type data_;

  BodyBind(const Tree& tree)
      : data(tree.bodies()) {
  }
  template <typename Range>
  BodyBind(const Tree& tree, const Range& range)
      : data(tree.body_permute(range.begin(), tree.body_begin()),
             tree.body_permute(range.begin(), tree.body_end())) {
  }

  T& operator[](const typename Tree::body_type& body) {
    return data[body.index()];
  }
  const T& operator[](const typename Tree::body_type& body) const {
    return data[body.index()];
  }

  iterator operator[](const typename Tree::body_iterator& bi) {
    return data.begin() + (*bi).index();
  }
  const_iterator operator[](const typename Tree::body_iterator& bi) const {
    return data.begin() + (*bi).index();
  }

  std::pair<iterator,iterator> operator[](const typename Tree::box_type& box) {
    return {(*this)[box.body_begin()], (*this)[box.body_end()]};
  }
};
