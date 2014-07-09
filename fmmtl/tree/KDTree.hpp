/* * * * * * * *
 * @file KDTree
 * @brief General class representing a KD-Tree
 **/

#include <vector>
#include <algorithm>

#include <iostream>
#include <iomanip>

#include <boost/range.hpp>
#include <boost/iterator/iterator_adaptor.hpp>

#include "fmmtl/numeric/Vec.hpp"

namespace fmmtl {
  
  template <unsigned DIM, template POINT>
  class KDTree {

    // PUBLIC PREDECLARATIONS
    public:
      
      // How to get point_type from PointIter?
      // XXX How to template on this without knowing type?
      typedef POINT point_type;

      // The type of this tree
      typedef KDTree tree_type;

      struct Box;
      typedef Box box_type;

      struct box_iterator;

    // PRIVATE METHODS
    private:

      /** @brief inserts all the points from p_first to p_last
       *    into the KDTree.
       *    Continuously partitions based on a Functor comp that splits pivots the 
       *    midpt.
       *
       * @param[in]: p_first, PointIter at the start of range
       * @param[in]: p_last, PointIter at the end of range
       * @param[in]: comp, a Functor that takes in a midpt value and returns a comparison
       *
       **/

       //XXX Should be done recursively
      template <typename PointIter, typename Comparator>
      void insert(PointIter p_first, PointIter p_last, Comparator comp, int idx = 0, int level = 0) {
       
        // maybe can use this to track leafs?
        if (p_first == p_last) {
          return;
        }

        // get distance between the iterators (# pts)
        int length = std::distance(pf, pl);

        // select the axis based on level so that the axis cycles 
        // through all valid dimensions
        int axis = level % DIM;
  
        // Get the pivot midpt
        PointIter p_midpt = std::advance(pf, length / 2);
        

        // Use std::nth element to rearrange elements in the range [first, last), in such
        // a way that the element in the 'nth' positon is the element that would be in that
        // position in a sorted sequence.  
        // No other specific order, except that none of the elements preceding 'nth' are greater 
        // than it, and none of the elements following it are less
        std::nth_element(pf, pm, pl, comp(axis));

        // Get the new pf and new pl since the old ones might be invalidated after
        // calling nth_element, since it could have been shuffled
        // XXX: possible off by 1 error??
        Point pf = std::advance(p_midpt, -1 * (length / 2));
        Point pl = std::advance(p_midpt, length / 2);
        
        // XXX: IDEA: since each box relation is binary, then given each box idx X, then its child
        // will be 2X and 2X + 1, and its parent will be X << 1
        // create a new box_data and track it in box_data_
        // XXX: How to track leafs
        box_data_.push(box_data(level, idx / 2, std::distance(pf, p_midpt), 
                  std::distance(p_midpt, pl), *p_midpt)

        insert(p_first, p_midpt, comp, idx * 2, level + 1);
        insert(p_midpt, p_last, comp, idx * 2 + 1, level + 1);
      }

    // PUBLIC METHODS
    public:

      // Public constructor for an invalid KDTree
      template <typename PointIter>
      KDTree() {
        // invalid tree  
      }

      // Destructor
      ~KDTree();

      // Public constructor for a KDTree given a first and last iter through points
      // XXX: Interface conflict: level?
      //    - I'm using it because I am creating trees recursively
      // @param[in]: first, PointIter to start of range
      // @param[in]: last, PointIter to end of range
      // @param[in]: comp, Functor that takes in level (int) and returns a unary comparator
      //    that compares based on the axis
      // @param[in]: level, optional that gives level of tree
      // 
      // @pre: std::distance(first, last) >= 0
      // @pre: for all level, comp(level) provides a valid comparator
      template <typename PointIter, typename Comparator>
      KDTree(PointIter first, PointIter last, Comparator comp, int level = 0) {
        
        // insert all the values from first to last, then partition
        insert(first, last, level, n_crit);
      }

    // PRIVATE STRUCTS AND FIELDS
    private:

      typedef <template PointIter>
      struct KDNode {
    
          PointIter position_;
          KDTree left_;
          KDTree right_;

          // Public Constructor
          KDNode (PointIter position, KDTree left, KDTree right) :
            position_(position), left_(left), right_(right)
          {}
      }

      struct box_data {
       
        // index into our PointIter vector
        int idx;

        // precomputed center
        point_type center_;

        box_data (int idx, tree_type* tree, point_type center)
          : idx_(idx), tree_(tree), center_(center) {

        }

        // Accessors 
        unsigned parent() const {
          return idx_ << 1;
        }

        // XXX: not sure how to handle is leaf, maybe just set??
        bool is_leaf() {
          return level() == tree_.max_depth();
        }

        unsigned level() {
          // log_2(x) == log_10(x) / log_10(2)
          return log10(idx_ + 1) / log10(2);
        }

      };

      struct Box {
        typedef typename tree_type::box_iterator box_iterator;
      
        Box()
          : idx_(-1), tree_(nullptr) {

        }

        unsigned index() const {
          return idx_;
        }

        unsigned level() const {
          return data().leve();
        }

        // point_type extents() const {}

        // double volume() const {}

        // double radius() const {}

        unsigned num_children() const {
          return std::distance(child_begin(), child_end());
        }

        // unsigned num_bodies() const {}

        bool is_leaf() const {
          return data().is_leaf();
        }

        point_type center() const {
          return data().center_;
        }

        Box parent() const {
          return Box(data().parent(), tree_);
        }

        // XXX: body_iterators
        box_iterator child_begin() const {
          FMMTL_ASSERT(!is_leaf());
          return box_iterator(data().begin_, tree_);
        }

        box_iterator child_end() const {
          FMMTL_ASSERT(!is_leaf());
          return box_iterator(data().end_, tree_);
        }

        bool operator==(const Box& b) const {
          
          // Why can't we just track trees??
          FMMTL_ASSERT(tree_ == b.tree_);
          return this->index() == b.index();
        }
    
        bool operator<(const Box& b) const {
          FMMTL_ASSERT(tree_ == b.tree_);
          return this->index() < b.index();
        }

      private:
        int idx_; // XXX: Cris uses unsigned, maybe we stick with that?
        tree_type* tree_;
        
        // Private Constructor for friends-only
        Box (int idx, const tree_type* tree)
          : idx_(idx), tree_(const_cast<tree_type*>(tree)) {
        }

        inline box_data& data() const {
          return tree_->box_data_[idx_];
        }
        
        friend class KDTree;

      };

      /** @struct KDTree::box_iterator 
       * @brief: box_iterator to iterate through boxes
       */
      struct box_iterator
        : public iterator_adaptor<box_iterator,
                                  unsigned
                                  Box,
                                  std::random_access_iterator_tag,
                                  Box,
                                  unsigned>
      {
        inline box_iterator()
          : box_iterator::iterator_adaptor(0), tree_(nullptr) {
        }

        private:
          const tree_type* tree_;
          friend class KDTree;
          
          // box_iterator Constructors
          // XXX: Implement proxy design
          inline box_iterator(unsigned idx, const tree_type* tree)
            : box_iterator::iterator_adaptor(idx), tree_(tree) {
          }

          inline box_iterator(const Box& b)
            : box_iterator::iterator_adaptor(b.index()), tree_(b.tree_) {
          }

          friend class boost::iterator_core_access;

          // iterator_adaptor override iterator traits
          inline Box dereference() const {
            return Box(index(), tree_);
          }

          inline unsigned index() const {
            return this->base_reference();
          }
      };

      point_type center() const {
        return *position_;
      }

      inline unsigned size() const {
        return permute_.size();
      }

      inline unsigned boxes() const {
        return box_data_.size();
      }

      inline unsigned levels() const {
        return level_offset_.size() - 1;
      }

      inline bool contains (const box_type& box) const {
        return this == box.tree_;
      }

      box_type root() const {
        return Box(0, this);
      }

      box_type box (const int index) const {
        FMMTL_ASSERT(idx < box_data_.size());
        return Box(idx, this);
      }

      box_iterator box_begin() const {
        return box_iterator(0, this);
      }

      box_iterator box_end() const {
        return box_iterator(boxes, this);
      }

      box iterator box_begin (int L) const {
        FMMTL_ASSERT(L < levels());
        return box_iterator (level_offset_[L], this);
      }

      box_iterator box_end (int L) const {
        FMMTL_ASSERT( L < levels());
        return box_iterator(level_offset_[l+1], this);
      }

      // PRIVATE FIELDS

      // Permutation: permute_[i] is the current idx of originally ith point
      std::vector<unsigned> permute_;
      // level_offset_[i] and level_offset_[i+1] is the start and end of a level[i]
      std::vector<unsigned> level_offset_;
      // Vector of data describing a box
      std::vector<box_data> box_data_;
      
      
      KDNode node_;
      int level_;


  }
