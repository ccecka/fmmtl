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
  
  template <typename DIM>
  class KDTree {
    public:

      // The spacial point type used for centers and extents
      typedef Vec<DIM, double> point_type;

      // The type of this tree
      typedef KDTree tree_type;

      // Predeclarations
      struct Body;
      typedef Body body_type;


    private:
      template <typename PointIter>
      void insert(PointIter p_first, PointIter p_last, unsigned NCRIT) {

        


      }


    public:

      // Public constructor for our KDTree
      template <typename PointIter>
      KDTree(PointIter first, PointIter last, unsigned n_crit = 256) {
        
        // insert all the values from first to last, then partition
        insert(first, last, n_crit);
      }


  }
