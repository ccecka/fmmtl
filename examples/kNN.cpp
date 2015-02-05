#include <iostream>
#include <vector>
#include <array>

#include "fmmtl/Direct.hpp"
#include "fmmtl/numeric/random.hpp"
#include "fmmtl/tree/KDTree.hpp"
#include "fmmtl/tree/TreeData.hpp"

#include "fmmtl/traversal/SingleTraversal.hpp"

template <typename T, std::size_t K, typename Compare = std::less<T> >
class ordered_vector {
  static_assert(K > 0, "ordered_vector must have K > 0");

  // Empty base class optimization for trivial comparators
  struct internal : public Compare {
    internal(const Compare& comp, unsigned size)
        : Compare(comp), size_(size) {}
    unsigned size_;
    std::array<T,K> data_;
  };
  internal int_;

  // Insert @a v in order into @a int_.data_
  void insert(unsigned i, const T& v) {
    for ( ; i > 0 && int_(v,int_.data_[i-1]); --i)
      int_.data_[i] = int_.data_[i-1];
    int_.data_[i] = v;
  }

 public:
  using value_type = T;
  using const_iterator = typename std::array<T,K>::const_iterator;

  // Default construction
  ordered_vector(const Compare& comp = Compare())
      : int_(comp, 0) {
  }

  unsigned size() const {
    return int_.size_;
  }
  const value_type& operator[](unsigned i) const {
    return int_.data_[i];
  }
  const_iterator begin() const {
    return int_.data_.begin();
  }
  const_iterator end() const {
    return int_.data_.begin() + int_.size_;
  }

  const T& back() const {
    return int_.data_[int_.size_-1];
  }

  ordered_vector& operator+=(const T& v) {
    if (size() < K)
      insert(int_.size_++, v);
    else if (int_(v,int_.data_[K-1]))
      insert(K-1, v);
    return *this;
  }

  operator std::vector<T>() const {
    return std::vector<T>(begin(), end());
  }

  bool operator==(const ordered_vector& v) const {
    return std::equal(begin(), end(), v.begin());
  }
};  // end ordered_vector


/** Print an ordered_vector to an output stream */
template <class T, std::size_t K, class C>
std::ostream& operator<<(std::ostream& s, const ordered_vector<T,K,C>& ov) {
  s << "(";
  auto first = ov.begin(), last = ov.end();
  if (first < last)
    s << *first;
  for (++first; first < last; ++first)
    s << ", " << *first;
  return s << ")";
}


/** kNN Kernel implementing
 * r_i += K(t_i, s_j) c_j
 * where K(t_i, s_j) = ||t_i - s_j||^2
 * r_i is a sorted vector of the K smallest distance-idx pairs seen
 * and c_j = j
 */
template <std::size_t K>
struct kNN {
  typedef Vec<3,double> source_type;
  typedef Vec<3,double> target_type;
  typedef unsigned      charge_type;

  struct dist_idx_pair {
    double distance_sq;
    unsigned index;
    bool operator<(const dist_idx_pair& other) const {
      return distance_sq < other.distance_sq;
    }
    bool operator==(const dist_idx_pair& other) const {
      return distance_sq == other.distance_sq && index == other.index;
    }
    friend std::ostream& operator<<(std::ostream& s, const dist_idx_pair& dip) {
      return s << "(" << dip.distance_sq << ", " << dip.index << ")";
    }
  };
  /** A kernel_value_type is the result of ||t_i - s_j||, a double.
   * Multiplication with a charge (unsigned),
   * is a concatenation into dist_idx_pair
   */
  struct kernel_value_type {
    double distance_sq;
    dist_idx_pair operator*(const charge_type& c) const {
      return {distance_sq, c};
    }
  };

  typedef ordered_vector<dist_idx_pair,K> result_type;

  kernel_value_type operator()(const target_type& t,
                               const source_type& s) const {
    return {norm_2_sq(t-s)};
  }
  kernel_value_type transpose(const kernel_value_type& kvt) const {
    return kvt;
  }
};



int main(int argc, char** argv) {
  int N = 1000;
  int M = 1000;
  bool checkErrors = true;
  Clock timer;

  // Parse custom command line args
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i],"-N") == 0) {
      N = atoi(argv[++i]);
    } else if (strcmp(argv[i],"-M") == 0) {
      M = atoi(argv[++i]);
    } else if (strcmp(argv[i],"-nocheck") == 0) {
      checkErrors = false;
    }
  }

  // Define the kernel types
  using Kernel = kNN<5>;
  using source_type = typename Kernel::source_type;
  using charge_type = typename Kernel::charge_type;
  using target_type = typename Kernel::target_type;
  using result_type = typename Kernel::result_type;

  // Define the tree types
  constexpr unsigned DS = fmmtl::dimension<source_type>::value;
  constexpr unsigned DT = fmmtl::dimension<target_type>::value;
  using SourceTree = fmmtl::KDTree<DS>;
  using TargetTree = fmmtl::KDTree<DT>;
  using source_box_type  = typename SourceTree::box_type;
  using source_body_type = typename SourceTree::body_type;
  using target_box_type  = typename TargetTree::box_type;
  using target_body_type = typename TargetTree::body_type;

  // Construct the kernel
  Kernel K;
  //std::cout << KernelTraits<Kernel>() << std::endl;

  // Construct the kernel data
  std::vector<source_type> sources = fmmtl::random_n(M);
  std::vector<charge_type> charges(M);

  std::vector<target_type> targets = fmmtl::random_n(N);
  std::vector<result_type> results(N);

  // Charges are the indices of the original sources
  std::iota(charges.begin(), charges.end(), 0);

  timer.start();
  // Construct the source tree
  SourceTree source_tree(sources);

  // Permute the sources and charges to the source_tree
  auto p_sources = make_body_binding(source_tree, sources);
  auto p_charges = make_body_binding(source_tree, charges);

  double construct_time = timer.seconds();
  std::cout << "Construct: " << construct_time << std::endl;

  //
  // Rules -- kNN Single-Tree
  //

  // Precompute the bounding box of each box
  using bounding_box_type = fmmtl::BoundingBox<typename SourceTree::point_type>;
  auto box_bb = make_box_binding<bounding_box_type>(source_tree);
  for (source_box_type b : boxes(source_tree))
    box_bb[b] = bounding_box_type(b.center() - b.extents()/2,
                                  b.center() + b.extents()/2);

  //
  // Traversal -- Single Tree
  //
  timer.start();

  // For each target
#pragma omp parallel for
  for (unsigned k = 0; k < targets.size(); ++k) {
    // Get the target and result
    const target_type& t = targets[k];
    result_type& r = results[k];

    // Associate each box of the source tree with a distance from target t
    auto hyper_rect = make_box_binding<double>(source_tree);
    double max_distance_sq = 1e200;

    // Define the rules of the traversal
    auto base = [&](const source_box_type& b) {
      // For all the sources/charges of this box
      auto ci = p_charges[b.body_begin()];
      for (const source_type& s : p_sources[b]) {
        r += K(t,s) * (*ci);
        ++ci;
      }
      max_distance_sq = r.back().distance_sq;
    };
    auto prune = [&](const source_box_type& b) {
      return hyper_rect[b] >= max_distance_sq;
    };
    auto visit_order = [&](const source_box_type& b) {
      // Make sure this is a binary tree
      static_assert(b.num_children() == 2, "Binary Tree Only For Now");

      source_box_type c_left  = *(  b.child_begin());
      source_box_type c_right = *(++b.child_begin());
      // Compute and store the distance to the target
      double left  = hyper_rect[c_left]  = norm_2_sq(box_bb[c_left], t);
      double right = hyper_rect[c_right] = norm_2_sq(box_bb[c_right], t);

      if (left < right) return std::array<source_box_type,2>{c_left, c_right};
      else              return std::array<source_box_type,2>{c_right, c_left};
    };

    // Traverse the source tree
    fmmtl::traverse(source_tree.root(), prune, base, visit_order);
  }

  double traverse_time = timer.seconds();
  std::cout << "Traverse: " << traverse_time << std::endl;

  //
  // Complete
  //

  // Check the result
  if (checkErrors) {
    std::cout << "Computing direct..." << std::endl;

    std::vector<result_type> exact(N);

    // Compute the result with a direct matrix-vector multiplication
    timer.start();
    fmmtl::direct(K, sources, charges, targets, exact);
    double direct_time = timer.seconds();
    std::cout << "Direct: " << direct_time << std::endl;

    int wrong_results = 0;
    for (unsigned k = 0; k < results.size(); ++k) {
      if (!(exact[k] == results[k])) {
        std::cout << "[" << std::setw(log10(M)+1) << k << "]"
                  << " Exact: " << exact[k]
                  << ", Tree: " << results[k] << std::endl;
        ++wrong_results;
      }
    }
    std::cout << "Wrong counts: " << wrong_results << " of " << M << std::endl;
  }

  return 0;
}
