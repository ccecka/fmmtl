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
  
  class KDTree {
    public:

      // The type of this tree
      typedef KDTree tree_type;

      // Predeclarations
      struct KDNode;
      typedef KDNode node_type;

    private:
      template <typename PointIter, typename Comparator>
      KDNode insert(PointIter p_first, PointIter p_last, Comparator comp, int depth, unsigned NCRIT) {
        
        // get distance between the iterators (# pts)
        int length = std::distance(p_first, p_last);

        // select the axis based on depth so that the axis cycles 
        // through all valid values
        int axis = depth % length;

        // Get the appropriate comparator using functor, and then
        // sorts iterators based on the axis
        std::sort(p_first, p_last, comp(axis));

        // XXX: How to get new start after sorting??
        PointIter p_midpt = std::advance(p_first, length / 2);

        depth_ = depth;
        node_ = KDNode (p_midpt, KDTree(p_first, p_midpt, comp, depth + 1, NCRIT), KDTree(std::advance(p_midpt, 1), p_last, comp, depth + 1, NCRIT));

      }


    public:

      // Public constructor for an invalid KDTree
      template <typename PointIter>
      KDTree() {
        // invalid tree  
      }

      // Destructor
      ~KDTree();

      // Public constructor for a KDTree given a first and last iter through points
      // XXX: Interface conflict: depth?
      //    - I'm using it because I am creating trees recursively
      // @param[in]: first, PointIter to start of range
      // @param[in]: last, PointIter to end of range
      // @param[in]: comp, Functor that takes in depth (int) and returns a unary comparator
      //    that compares based on the axis
      // @param[in]: depth, optional that gives depth of tree
      // 
      // @pre: std::distance(first, last) >= 0
      // @pre: for all depth, comp(depth) provides a valid comparator
      template <typename PointIter, typename Comparator>
      KDTree(PointIter first, PointIter last, Comparator comp, int depth = 0) {
        
        // insert all the values from first to last, then partition
        insert(first, last, depth, n_crit);
      }

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

      // PRIVATE FIELDS
      KDNode node_;
      int depth_;

  }
