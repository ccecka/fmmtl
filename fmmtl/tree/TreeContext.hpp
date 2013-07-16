#pragma once

template <class Tree>
struct TreeTraits {
  typedef typename Tree::box_type      box_type;
  typedef typename Tree::box_iterator  box_iterator;
  typedef typename Tree::body_type     body_type;
  typedef typename Tree::body_iterator body_iterator;
};

template <class DerivedTree>
struct Tree;

template <class DerivedTree>
struct Box;

template <class DerivedTree>
struct Body;



template <class P>
struct Tree<Octree<P>> {
  typedef typename TreeTraits<Tree>::box_type      box_type;
  typedef typename TreeTraits<Tree>::box_iterator  box_iterator;
  typedef typename TreeTraits<Tree>::body_type     body_type;
  typedef typename TreeTraits<Tree>::body_iterator body_iterator;

  // Tree interfacing with Octree
};


template <class P>
struct Body<Octree<P>> {
  // Body interfacing with Octree::body_type
};

template <class P>
struct Box<Octree<P>> {
  // Box interfacing with Octree::box_type
  typename DerivedTree::point_type extents() const {
    return static_cast<const Derived*>(this)->extents();
  }
};



template <class Tree, class Data>
class BoxData;

template <class Tree, class Data>
class BodyData;

template <class Tree, class Data>
class BodyAssoc;



#include "Octree.hpp"

template <class P,
          class Data>
class BoxData<Octree<P>, Data> {
  //! Data corresponding to box indices in a tree
  typedef std::vector<Data> data_container;
  data_container data_;

 public:
  inline Data& operator[](const typename Tree::box_type& box) {
    return data_[box.index()];
  }
  inline const Data& operator[](const typename Tree::box_type& box) const {
    return data_[box.index()];
  }
};

// Associate a range of data with the bodies in the tree.
// The data is assumed to be sen in the same order as the tree was initially
// constructed.
template <class P,
          class Data>
class BodyData<Octree<P>, Data> {
  //! Data corresponding to body indices in a tree
  typedef std::vector<Data> data_container;
  data_container data_;

 public:
  template <class IT>
  BodyAssoc(IT first, IT last)
      :

  inline Data& operator[](const typename Tree::body_type& body) {
    return data_[body.index()];
  }
  inline const Data& operator[](const typename Tree::body_type& body) const {
    return data_[body.index()];
  }

  inline ITER begin() {
    return data_.begin();
  }
  inline PITER begin(const typename Tree::body_iterator& bi) {

  }
  inline PITER begin(const typename Tree::box_type& box) {
    return begin(box.body_begin());
  }
};
