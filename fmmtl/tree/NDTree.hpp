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
#include "fmmtl/tree/MortonCoder.hpp"

namespace fmmtl {
using boost::has_range_iterator;
using boost::iterator_adaptor;

/** In-place bucket sort using counting sort
 *
 * @param[in] first,last Iterator pair to the sequence to be bucketed
 * @param[in] num_buckets The number of buckets to be used
 * @param[in] map Mapping of elements in [first,last) to [0,num_buckets)
 * @returns vector of iterators:
 *          [result[i],result[i+1]) is the range of the ith bucket.
 *
 * @pre For all i in [first,last), 0 <= map(*i) < num_buckets.
 * @post For all i and j such that first <= i < j < last,
 *       then 0 <= map(*i) <= map(*j) < num_buckets.
 * @tparam Iterator models a random-access mutable iterator
 */
template <typename Iterator, typename BucketMap>
std::vector<Iterator> bucket_sort(Iterator first, Iterator last,
                                  unsigned num_buckets, BucketMap map) {
  typedef typename std::iterator_traits<Iterator>::value_type value_type;

  std::vector<Iterator> start(num_buckets+1, first);

  for ( ; first != last; ++first)
    ++start[1 + map(*first)];

  for (unsigned k = 2; k <= num_buckets; ++k)
    start[k] += (start[k-1] - start[0]);

  std::vector<Iterator> end = start;

  // For each bin
  for (unsigned curr_bin = 0; curr_bin < num_buckets; ++curr_bin) {
    first = start[curr_bin];
    last  = end[1+curr_bin];

    while (first != last) {
      value_type& v = *first++;

      // Swap until an element comes back to this bin
      unsigned bin;
      while (curr_bin != (bin = map(v)))
        std::swap(v, *start[bin]++);
    }
  }

  return end;
}


//! Class for tree structure
template <unsigned DIM>
class NDTree {
 public:
  //! The spacial point type used for centers and extents
  typedef Vec<DIM,double> point_type;

  //! Each box has 2^DIM children
  static constexpr unsigned max_children = 1 << DIM;

  //! The type of this tree
  typedef NDTree<DIM> tree_type;

  // Predeclarations
  struct Body;
  typedef Body body_type;
  struct Box;
  typedef Box box_type;
  struct body_iterator;
  struct box_iterator;

 private:
  // Morton coder to use for the points
  MortonCoder<DIM> coder_;
  // Code type
  typedef typename MortonCoder<DIM>::code_type code_type;

  struct box_data {
    // The level of the box and the leaf_bit
    unsigned level_;
    // The index of the parent of this box
    unsigned parent_;
    // These can be either point offsets or box offsets depending on is_leaf
    unsigned begin_;
    unsigned end_;

    // Precomputed center
    point_type center_;

    static constexpr unsigned leaf_bit = (1 << (8*sizeof(level_)-1));

    box_data(unsigned level, unsigned parent,
             unsigned begin, unsigned end,
             const point_type& center)
        : level_(level), parent_(parent),
          begin_(begin), end_(end),
          center_(center) {
    }
    unsigned parent() const {
      return parent_;
    }
    void set_leaf() {
      level_ |= leaf_bit;
    }
    void set_non_leaf() {
      level_ &= ~leaf_bit;
    }
    bool is_leaf() const {
      return level_ & leaf_bit;
    }
    unsigned level() const {
      return level_ & ~leaf_bit;
    }
  };

  // Morton coded objects this Tree holds.
  //std::vector<point_type> point_;

  // Morton code for each point
  std::vector<code_type> mc_;
  // Permutation: permute_[i] is the current idx of originally ith point
  std::vector<unsigned> permute_;
  // level_offset_[i] and level_offset_[i+1] is the start and end of level i
  std::vector<unsigned> level_offset_;
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
    code_type morton_index() const {
      return tree_->mc_[idx_];
    }
   private:
    unsigned idx_;
    tree_type* tree_;
    Body(unsigned idx, const tree_type* tree)
        : idx_(idx), tree_(const_cast<tree_type*>(tree)) {
      FMMTL_ASSERT(idx_ < tree_->size());
    }
    friend class NDTree;
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
      return data().level();
    }
    point_type extents() const {
      const BoundingBox<point_type> bb = tree_->coder_.bounding_box();
      return (bb.max() - bb.min()) / (1 << level());
    }
    double volume() const {
      point_type e = extents();
      return std::accumulate(e.begin(), e.end(), double(1), std::multiplies<double>());
    }
    double radius() const {
      return norm_2(extents()) / 2.0;
    }
    unsigned num_children() const {
      return std::distance(child_begin(), child_end());
    }
    unsigned num_bodies() const {
      return std::distance(body_begin(), body_end());
    }
    bool is_leaf() const {
      return data().is_leaf();
    }
    /** The center of this box */
    point_type center() const {
      return data().center_;
    }
    /** The parent box of this box */
    Box parent() const {
      return Box(data().parent(), tree_);
    }

    /** The begin iterator to the bodies contained in this box */
    body_iterator body_begin() const {
      if (is_leaf())
        return body_iterator(data().begin_, tree_);
      else
	      return child_begin()->body_begin();
    }
    /** The end iterator to the bodies contained in this box */
    body_iterator body_end() const {
      if (is_leaf())
        return body_iterator(data().end_, tree_);
      else
	      return (--child_end())->body_end();
    }

    /** The begin iterator to the child boxes contained in this box */
    box_iterator child_begin() const {
      FMMTL_ASSERT(!is_leaf());
      return box_iterator(data().begin_, tree_);
    }
    /** The end iterator to the child boxes contained in this box */
    box_iterator child_end() const {
      FMMTL_ASSERT(!is_leaf());
      return box_iterator(data().end_, tree_);
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
    friend class NDTree;
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
    inline unsigned index() const {
      return this->base_reference();
    }
   private:
    const tree_type* tree_;
    friend class NDTree;
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
    inline unsigned index() const {
      return this->base_reference();
    }
   private:
    const tree_type* tree_;
    friend class NDTree;
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
  };

  /** Construct an tree encompassing a bounding box
   * and insert a range of points */
  template <typename Range>
  NDTree(const Range& rng, unsigned n_crit = 256,
         typename std::enable_if<has_range_iterator<Range>::value>::type* = 0)
      : NDTree(rng.begin(), rng.end(), n_crit) {
  }

  /** Construct an tree encompassing a bounding box
   * and insert a range of points */
  template <typename PointIter>
  NDTree(PointIter first, PointIter last, unsigned n_crit = 256)
      : coder_(get_boundingbox(first, last)) {
    insert(first, last, n_crit);
  }

  /** Return the Bounding Box that this NDTree encompasses */
  BoundingBox<point_type> bounding_box() const {
    return coder_.bounding_box();
  }

  /** Return the center of this NDTree */
  point_type center() const {
    return coder_.center();
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

  /** The number of boxes contained in level L of this tree */
  inline unsigned boxes(unsigned L) const {
    return level_offset_[L+1] - level_offset_[L];
  }

  /** The maximum level of any box in this tree */
  inline unsigned levels() const {
    return level_offset_.size() - 1;
  }
  /** The maximum possible level of any box in this tree */
  inline static unsigned max_level() {
    return MortonCoder<DIM>::levels() - 1;
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
    FMMTL_ASSERT(L < levels());
    return box_iterator(level_offset_[L], this);
  }
  /** Return an iterator one past the last box at level L in this tree
   * @pre L < levels()
   */
  box_iterator box_end(unsigned L) const {
    FMMTL_ASSERT(L < levels());
    return box_iterator(level_offset_[L+1], this);
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
  body_permute(RandomAccessIter&& it, const body_iterator& bi) const {
    return boost::make_permutation_iterator(it, permute_.cbegin() + bi.index());
  }

  /** Write an NDTree to an output stream */
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
    FMMTL_LOG("Tree Insert");

    // Create a code-idx pair vector
    typedef std::pair<code_type, unsigned> code_pair;
    std::vector<code_pair> codes;
    codes.reserve(std::distance(p_first, p_last));
    //std::vector<point_type> points;

    unsigned idx = 0;
    for (PointIter pi = p_first; pi != p_last; ++pi, ++idx) {
      point_type p = *pi;
      FMMTL_ASSERT(coder_.bounding_box().contains(p));

      //points.push_back(p);
      codes.emplace_back(coder_.code(p), idx);
    }

    // Push the root box which contains all points
    box_data_.emplace_back(0, 0, 0, codes.size(), center());
    level_offset_.push_back(0);

    // For every box that is created
    for (unsigned k = 0; k != box_data_.size(); ++k) {
      // @pre box_data_[k] has not been designated a leaf yet
      // @pre box_data_[k].begin_ and end_ refer to body indices

      // If this box has few enough points, mark as leaf and go to next box
      if (box_data_[k].end_ - box_data_[k].begin_ <= NCRIT) {
        box_data_[k].set_leaf();
        continue;
      }
      // Else, split this box

      // If the children will start a new level, record it
      if (box_data_[k].level() + 1 > levels())
        level_offset_.push_back(box_data_.size());

      // Get the box data
      auto code_begin = codes.begin() + box_data_[k].begin_;
      auto code_end   = codes.begin() + box_data_[k].end_;
      const unsigned shift = DIM * (max_level() - box_data_[k].level());

      // Sort the points in this box into the "bucket" children
      auto off = bucket_sort(code_begin, code_end, max_children,
                             [=] (const code_pair& v)
                             { return (v.first >> shift) & (max_children-1); });

      // Record the child begin idx
      box_data_[k].begin_ = box_data_.size();

      // For each bucket
      for (unsigned c = 0; c < max_children; ++c) {
        // If this child contains points
        if (off[c+1] != off[c]) {
          // Add the child
          unsigned level = box_data_[k].level() + 1;
          box_data_.emplace_back(level,                      // Level
                                 k,                          // Parent idx
                                 off[c]   - codes.begin(),   // Body begin idx
                                 off[c+1] - codes.begin(),   // Body end idx
                                 get_center(off[c]->first, level));  // Center
        }
      }

      // Record the child end idx
      box_data_[k].end_ = box_data_.size();
    }

    // Record the end of the last level
    level_offset_.push_back(box_data_.size());

    // Allocate
    mc_.reserve(codes.size());
    permute_.reserve(codes.size());
    // Extract the code, permutation vector, and sorted point
    for (auto& c : codes) {
      mc_.push_back(c.first);
      permute_.push_back(c.second);
      //point_.push_back(points[permute_.back()]);
    }
  }

  /** Get the center of the box that
   * morton code @a c is contained in at level @level
   */
  point_type get_center(code_type c, unsigned level) {
    // Mask for boxes of this level
    code_type mask = code_type(1) << (DIM*(max_level() - level + 1));
    --mask;
    return coder_.center(c & ~mask /*cmin*/, c | mask /*cmax*/);
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
  NDTree(const NDTree&) {};
  void operator=(const NDTree&) {};
};

} // end namespace fmmtl
