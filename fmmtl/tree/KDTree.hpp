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
  
  template <unsigned DIM>
  class KDTree {

    // PUBLIC PREDECLARATIONS
    public:

      // The type of this tree
      typedef KDTree tree_type;

      // Predeclarations
      struct KDNode;
      typedef KDNode node_type;

      struct Box;
      typedef Box box_type;

      struct box_iterator;

    // PRIVATE METHODS
    private:
      template <typename PointIter, typename Comparator>
      KDNode insert(PointIter p_first, PointIter p_last, Comparator comp, int level, unsigned NCRIT) {
        
        // get distance between the iterators (# pts)
        int length = std::distance(p_first, p_last);

        // select the axis based on level so that the axis cycles 
        // through all valid dimensions
        int axis = level % DIM;

        // Get the appropriate comparator using functor, and then
        // sorts iterators based on the axis
        std::sort(p_first, p_last, comp(axis));

        // XXX: How to get new start after sorting??
        PointIter p_midpt = std::advance(p_first, length / 2);

        level_ = level;
        node_ = KDNode (p_midpt, KDTree(p_first, p_midpt, comp, level + 1, NCRIT), KDTree(std::advance(p_midpt, 1), p_last, comp, level + 1, NCRIT));

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

      // PRIVATE FIELDS
      KDNode node_;
      int level_;

  }
