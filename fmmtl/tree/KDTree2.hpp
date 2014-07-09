#pragma once
/** @file NDTree
 * @brief General class representing a {1D,2D,3D,4D}-Tree.
 */

#include <vector>
#include <algorithm>

#include <iostream>
#include <iomanip>

#include <boost/range.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/permutation_iterator.hpp>

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/tree/BoundingBox.hpp"

namespace fmmtl {
using boost::has_range_iterator;
using boost::iterator_adaptor;


//! Class for tree structure
template <unsigned DIM>
class KDTree {
 public:
  //! The spacial point type used for centers and extents
  typedef Vec<DIM,double> point_type;

  //! Each box has 2 children
  static constexpr unsigned max_children = 2;

  //! The type of this tree
  typedef KDTree<DIM> tree_type;

  // Predeclarations
  struct Body;
  typedef Body body_type;
  struct Box;
  typedef Box box_type;
  struct body_iterator;
  struct box_iterator;

 private:

  // XXX: Dummy for now
  struct box_data {
    // Range (into permute_) of bodies this box contains
    unsigned begin, end;

    // Precomputed center
    //point_type center_;

    void set_leaf() {}  // XXX?
  };

  // Permutation: permute_[i] is the current idx of originally ith point
  std::vector<unsigned> permute_;
  // Vector of data describing a box
  std::vector<box_data> box_data_;

 public:

  struct Body {
    /** Construct an invalid Body */
    Body()
        : idx_(0), tree_(nullptr) {
    }
    //const point_type& point() const {
    //  return tree_->point_[idx_];
    //}
    //! The original order this body was seen
    unsigned number() const {
      return tree_->permute_[idx_];
    }
    //! The current order of this body
    unsigned index() const {
      return idx_;
    }
   private:
    unsigned idx_;
    tree_type* tree_;
    Body(unsigned idx, const tree_type* tree)
        : idx_(idx), tree_(const_cast<tree_type*>(tree)) {
      FMMTL_ASSERT(idx_ < tree_->size());
    }
    friend class KDTree;
  };

  // A tree-aligned box
  struct Box {
    typedef typename tree_type::box_iterator  box_iterator;
    typedef typename tree_type::body_iterator body_iterator;

    /** Construct an invalid Box */
    Box()
        : idx_(0), tree_(nullptr) {
    }
    unsigned index() const {
      return idx_;
    }
    unsigned level() const {
      return 0;  // XXX
    }
    point_type extents() const {
      return point_type(); // XXX
    }
    double volume() const {
      return 0; // XXX
    }
    double radius() const {
      return 0; // XXX
    }
    unsigned num_children() const {
      return 2;
    }
    unsigned num_bodies() const {
      return std::distance(body_begin(), body_end());
    }
    bool is_leaf() const {
      return num_bodies() < 256; // XXX
    }
    /** The center of this box */
    point_type center() const {
      return point_type();  // XXX
    }
    /** The parent box of this box */
    Box parent() const {
      return Box((idx_-1)/2, tree_);
    }

    /** The begin iterator to the bodies contained in this box */
    body_iterator body_begin() const {
      return body_iterator(data().begin, tree_);
    }
    /** The end iterator to the bodies contained in this box */
    body_iterator body_end() const {
      return body_iterator(data().end, tree_);
    }

    /** The begin iterator to the child boxes contained in this box */
    box_iterator child_begin() const {
      return box_iterator(2*idx_+1, tree_);
    }
    /** The end iterator to the child boxes contained in this box */
    box_iterator child_end() const {
      FMMTL_ASSERT(!is_leaf());
      return box_iterator(2*idx_+3, tree_);
    }

    /** Comparison operators for std:: algorithms */
    bool operator==(const Box& b) const {
      FMMTL_ASSERT(tree_ == b.tree_);
      return this->index() == b.index();
    }
    bool operator<(const Box& b) const {
      FMMTL_ASSERT(tree_ == b.tree_);
      return this->index() < b.index();
    }

    /** Write a Box to an output stream */
    inline friend std::ostream& operator<<(std::ostream& s,
                                           const box_type& b) {
      unsigned num_bodies = b.num_bodies();
      unsigned first_body = b.body_begin()->index();
      unsigned last_body = first_body + num_bodies - 1;

      return s << "Box " << b.index()
               << " (L" << b.level() << ", P" << b.parent().index()
               << ", " << num_bodies << (num_bodies == 1 ? " body" : " bodies")
               << " " << first_body << "-" << last_body
               << "): " << b.center() << " - " << b.extents();
    }
   private:
    unsigned idx_;
    tree_type* tree_;
    Box(unsigned idx, const tree_type* tree)
        : idx_(idx), tree_(const_cast<tree_type*>(tree)) {
    }
    inline box_data& data() const {
      return tree_->box_data_[idx_];
    }
    friend class KDTree;
  };

  /** @struct Tree::box_iterator
   * @brief Random access iterator for Boxes in the tree
   */
  struct box_iterator
      : public iterator_adaptor<box_iterator,                    // Derived
                                unsigned,                        // BaseType
                                Box,                             // Value
                                std::random_access_iterator_tag, // IterCategory
                                Box,                             // Reference
                                unsigned>                        // DiffType
  {
    /** Construct an invalid box_iterator */
    inline box_iterator()
        : box_iterator::iterator_adaptor(0), tree_(nullptr) {
    }
   private:
    const tree_type* tree_;
    friend class KDTree;
    inline box_iterator(unsigned idx, const tree_type* tree)
        : box_iterator::iterator_adaptor(idx), tree_(tree) {
    }
    inline box_iterator(const Box& b)
        : box_iterator::iterator_adaptor(b.index()), tree_(b.tree_) {
    }
    friend class boost::iterator_core_access;
    inline Box dereference() const {
      return Box(index(), tree_);
    }
    inline unsigned index() const {
      return this->base_reference();
    }
  };

  /** @struct Tree::body_iterator
   * @brief Random access iterator class for Bodies in the tree
   */
  struct body_iterator
      : public iterator_adaptor<body_iterator,                   // Derived
                                unsigned,                        // BaseType
                                Body,                            // Value
                                std::random_access_iterator_tag, // IterCategory
                                Body,                            // Reference
                                unsigned>                        // DiffType
  {
    /* Construct an invalid body_iterator */
    inline body_iterator()
        : body_iterator::iterator_adaptor(0), tree_(nullptr) {
    }
   private:
    const tree_type* tree_;
    friend class KDTree;
    inline body_iterator(unsigned idx, const tree_type* tree)
        : body_iterator::iterator_adaptor(idx), tree_(tree) {
    }
    inline body_iterator(const Body& b)
        : body_iterator::iterator_adaptor(b.index()), tree_(b.tree_) {
    }
    friend class boost::iterator_core_access;
    inline Body dereference() const {
      return Body(index(), tree_);
    }
    inline unsigned index() const {
      return this->base_reference();
    }
  };

  /** Construct an tree encompassing a bounding box
   * and insert a range of points */
  template <typename Range>
  KDTree(const Range& rng, unsigned n_crit = 256,
         typename std::enable_if<has_range_iterator<Range>::value>::type* = 0)
      : KDTree(rng.begin(), rng.end(), n_crit) {
  }

  /** Construct an tree encompassing a bounding box
   * and insert a range of points */
  template <typename PointIter>
  KDTree(PointIter first, PointIter last, unsigned n_crit = 256) {
    insert(first, last, n_crit);
  }

  /** Return the Bounding Box that this KDTree encompasses */
  BoundingBox<point_type> bounding_box() const {
    return BoundingBox<point_type>();  // XXX
  }

  /** Return the center of this KDTree */
  point_type center() const {
    return point_type(); // XXX
  }

  /** The number of bodies contained in this tree */
  inline unsigned size() const {
    return permute_.size();
  }
  /** The number of bodies contained in this tree */
  inline unsigned bodies() const {
    return size();
  }

  /** The number of boxes contained in this tree */
  inline unsigned boxes() const {
    return box_data_.size();
  }

  /** The maximum level of any box in this tree */
  inline unsigned levels() const {
    return 0;  // XXX
  }
  /** The maximum possible level of any box in this tree */
  inline static unsigned max_level() {
    return 10; // XXX
  }

  /** Returns true if the box is contained in this tree, false otherwise */
  inline bool contains(const box_type& box) const {
    return this == box.tree_;
  }
  /** Returns true if the body is contained in this tree, false otherwise */
  inline bool contains(const body_type& body) const {
    return this == body.tree_;
  }

  /** Return the root box of this tree */
  box_type root() const {
    return Box(0, this);
  }
  /** Return a box given it's index */
  box_type box(const unsigned idx) const {
    FMMTL_ASSERT(idx < box_data_.size());
    return Box(idx, this);
  }
  /** Return an iterator to the first body in this tree */
  body_iterator body_begin() const {
    return body_iterator(0, this);
  }
  /** Return an iterator one past the last body in this tree */
  body_iterator body_end() const {
    return body_iterator(size(), this);
  }
  /** Return an iterator to the first box in this tree */
  box_iterator box_begin() const {
    return box_iterator(0, this);
  }
  /** Return an iterator one past the last box in this tree */
  box_iterator box_end() const {
    return box_iterator(boxes(), this);
  }
  /** Return an iterator to the first box at level L in this tree
   * @pre L < levels()
   */
  box_iterator box_begin(unsigned L) const {
    return box_iterator((1 << L) - 1, this);
  }
  /** Return an iterator one past the last box at level L in this tree
   * @pre L < levels()
   */
  box_iterator box_end(unsigned L) const {
    return box_begin(L+1);
  }

  template <typename RandomAccessIter>
  struct body_permuted_iterator {
    typedef typename std::vector<unsigned>::const_iterator permute_iter;
    typedef boost::permutation_iterator<RandomAccessIter, permute_iter> type;
  };

  /** Tranform (permute) an iterator so its traversal follows the same order as
   * the bodies contained in this tree
   */
  template <typename RandomAccessIter>
  typename body_permuted_iterator<RandomAccessIter>::type
  body_permute(RandomAccessIter& it, const body_iterator& bi) const {
    return boost::make_permutation_iterator(it, permute_.cbegin() + bi.index());
  }

  /** Write an KDTree to an output stream */
  inline friend std::ostream& operator<<(std::ostream& s,
                                         const tree_type& t) {
    struct {
      inline std::ostream& print(std::ostream& s,
                                 const box_type& box) {
        s << std::string(2*box.level(), ' ') << box;
        if (!box.is_leaf())
          for (auto ci = box.child_begin(); ci != box.child_end(); ++ci)
            print(s << "\n", *ci);
        return s;
      }
    } recursive_box;

    return recursive_box.print(s, t.root());
  }

 private:
  //! TODO: Make dynamic and public?
  //! Uses incremental bucket sorting
  template <typename PointIter>
  void insert(PointIter p_first, PointIter p_last, unsigned NCRIT) {
    FMMTL_LOG("KDTree Insert");

    // Create a point-idx pair vector
    typedef typename std::iterator_traits<PointIter>::value_type point_i_type;
    typedef std::pair<point_i_type, unsigned> point_t;

    std::vector<point_t> point;
    point.reserve(std::distance(p_first, p_last));
    unsigned idx = 0;
    for (PointIter pi = p_first; pi != p_last; ++pi, ++idx) {
      point.emplace_back(*pi, idx);
    }

    // Push the root box which contains all points
    box_data_.push_back(box_data{0, unsigned(point.size())});

    // For every box that is created
    for (unsigned k = 0; k != box_data_.size(); ++k) {
      // @pre box_data_[k] has not been designated a leaf yet

      box_type box(k, this);

      if (box.num_bodies() <= NCRIT) {
        box_data_[k].set_leaf();
        continue;
      }

      // Else, split this box

      // Make a comparator for some dimension
      const unsigned dim = 0;   // Cycle or whatev
      auto comp = [](const point_t& a, const point_t& b) {
        return a.first[dim] < b.first[dim];
      };

      // Partition the points
      unsigned mid = box_data_[k].begin + (box_data_[k].end - box_data_[k].begin)/2;
      auto p_begin = point.begin() + box_data_[k].begin;
      auto p_mid   = point.begin() + mid;
      auto p_end   = point.begin() + box_data_[k].end;
      std::nth_element(p_begin, p_mid, p_end, comp);

      // Record the child boxes
      box_data_.push_back(box_data{box_data_[k].begin, mid});
      box_data_.push_back(box_data{mid, box_data_[k].end});
    }

    // Allocate
    permute_.reserve(point.size());
    // Extract the permutation idx
    for (auto&& p : point) {
      permute_.push_back(p.second);
    }
  }

  template <typename PointIter>
  BoundingBox<point_type> get_boundingbox(PointIter first, PointIter last) {
    // Construct a bounding box
    BoundingBox<point_type> bb(first, last);
    // Determine the size of the maximum side
    point_type extents =  bb.max() - bb.min();
    double max_side = *std::max_element(extents.begin(), extents.end());
    double radius = (1.0+1e-6) * max_side / 2.0;
    point_type center =  (bb.max() + bb.min()) / 2;
    // Make it square and add some wiggle room   TODO: Generalize on square
    return BoundingBox<point_type>(center - radius, center + radius);
  }

  // Just making sure for now
  KDTree(const KDTree&) {};
  void operator=(const KDTree&) {};
};

} // end namespace fmmtl
