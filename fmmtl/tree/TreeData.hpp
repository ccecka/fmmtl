#pragma once

#include <iterator>
#include <vector>

#include <thrust/iterator/permutation_iterator.h>


namespace fmmtl {

/** Maps boxes and box iterators to data and data iterators
 */
template <typename T, typename Tree>
struct BoxBind {
  typedef std::vector<T> container_type;
  typedef typename container_type::iterator iterator;
  typedef typename container_type::const_iterator const_iterator;

  container_type data;

  BoxBind() {
  }
  BoxBind(const Tree& tree)
      : data(tree.boxes()) {
  }

  bool empty() {
    return data.empty();
  }

  T& operator[](const typename Tree::box_type& box) {
    return data[box.index()];
  }
  const T& operator[](const typename Tree::box_type& box) const {
    return data[box.index()];
  }

  iterator operator[](const typename Tree::box_iterator& bi) {
    return std::begin(data) + (*bi).index();
  }
  const_iterator operator[](const typename Tree::box_iterator& bi) const {
    return std::begin(data) + (*bi).index();
  }

  iterator       begin()       { return std::begin(data); }
  iterator       end()         { return std::end(data); }
  const_iterator begin() const { return std::begin(data); }
  const_iterator end()   const { return std::end(data); }
};

template <typename T, typename Tree>
BoxBind<T,Tree> make_box_binding(const Tree& tree) {
  return {tree};
}


/** Permute iterators into the same order as in a tree
 */
template <typename Tree, typename RandomAccessIter>
using body_permute_iterator = thrust::permutation_iterator<RandomAccessIter,
                                                           typename Tree::permute_iterator>;

/** Tranform (permute) an iterator so its traversal follows the same order as
 * the bodies contained in a tree
 */
template <typename Tree, typename RandomAccessIter>
body_permute_iterator<Tree, RandomAccessIter>
permute_begin(const Tree& tree, RandomAccessIter iter) {
  return thrust::make_permutation_iterator(iter, tree.permute_begin());
}

template <typename Tree, typename RandomAccessIter>
body_permute_iterator<Tree, RandomAccessIter>
permute_end(const Tree& tree, RandomAccessIter iter) {
  return thrust::make_permutation_iterator(iter, tree.permute_end());
}

/** Tranform (permute) an iterator so its traversal follows the same order as
 * the bodies contained in a tree, starting at a particular body.
 *
 * @post
 *  permute_iter(tree, tree.body_begin(), iter) == permute_begin(tree, iter)
 *  permute_iter(tree, tree.body_end(),   iter) == permute_end(  tree, iter)
 */
template <typename Tree, typename RandomAccessIter>
body_permute_iterator<Tree, RandomAccessIter>
permute_iter(const Tree& tree,
             typename Tree::body_iterator bit,
             RandomAccessIter iter) {
  const auto i = std::distance(tree.body_begin(), bit);
  return thrust::make_permutation_iterator(iter, tree.permute_begin() + i);
}


/** Maps bodies and body iterators to data and data iterators
 */
template <typename T, typename Tree>
struct BodyBind {
  typedef std::vector<T> container_type;
  typedef typename container_type::iterator iterator;
  typedef typename container_type::const_iterator const_iterator;

  //const Tree& tree;
  container_type data;

  struct BodyDataRange {
    iterator b, e;
    iterator       begin()       { return b; }
    const_iterator begin() const { return b; }
    iterator       end()         { return e; }
    const_iterator end()   const { return e; }
    std::size_t    size()  const { return e - b; }
  };

  // Construct without permuting data or initialization data
  BodyBind(std::size_t n)
      : data(n) {
  }
  // Construct by permuting some initialization data based on the tree
  template <typename Iterator>
  BodyBind(const Tree& tree, Iterator data_it)
      : data(permute_begin(tree, data_it), permute_end(tree, data_it)) {
  }

  iterator       begin()       { return std::begin(data); }
  const_iterator begin() const { return std::begin(data); }
  iterator       end()         { return std::end(data); }
  const_iterator end()   const { return std::end(data); }

  T& operator[](const typename Tree::body_type& body) {
    return data[body.index()];
  }
  const T& operator[](const typename Tree::body_type& body) const {
    return data[body.index()];
  }
  T& operator[](typename Tree::size_type i) {
    return data[i];
  }
  const T& operator[](typename Tree::size_type i) const {
    return data[i];
  }

  iterator operator[](const typename Tree::body_iterator& bi) {
    return std::begin(data) + bi.index();
  }
  const_iterator operator[](const typename Tree::body_iterator& bi) const {
    return std::begin(data) + bi.index();
  }

  BodyDataRange operator[](const typename Tree::box_type& box) {
    return {operator[](box.body_begin()), operator[](box.body_end())};
  }
};

// Specialization for iterators
template <typename Tree, typename Iterator>
BodyBind<typename std::iterator_traits<Iterator>::value_type, Tree>
make_body_binding(const Tree& tree, Iterator range) {
  return {tree, range};
}

// Any Range with std::begin
template <typename Tree, typename Range>
auto
make_body_binding(const Tree& tree, const Range& range)
    -> decltype(make_body_binding(tree, std::begin(range))) {
  return make_body_binding(tree, std::begin(range));
}

// Any Type without initial data (avoid permutation, just allocate)
template <typename Type, typename Tree>
BodyBind<Type, Tree>
make_body_binding(const Tree& tree) {
  return {tree.bodies()};
}

} // end namespace fmmtl
