#pragma once

template <typename Tree>
struct TreeTraits {
  typedef Tree                              tree_type;
  typedef typename tree_type::box_type      box_type;
  typedef typename tree_type::body_type     body_type;
  typedef typename tree_type::box_iterator  box_iterator;
  typedef typename tree_type::body_iterator body_iterator;
};

template <typename SourceTree,
          typename TargetTree = SourceTree>
struct TreePairTraits {
  typedef TreeTraits<SourceTree> source_traits;
  typedef TreeTraits<TargetTree> target_traits;

  //! Source tree type
  typedef typename source_traits::tree_type        source_tree_type;
  //! Target tree type
  typedef typename target_traits::tree_type        target_tree_type;

  //! Source tree box type
  typedef typename source_traits::box_type      source_box_type;
  //! Source tree box iterator
  typedef typename source_traits::box_iterator  source_box_iterator;
  //! Source tree body type
  typedef typename source_traits::body_type     source_body_type;
  //! Source tree body iterator
  typedef typename source_traits::body_iterator source_body_iterator;

  //! Target tree box type
  typedef typename target_traits::box_type      target_box_type;
  //! Target tree box iterator
  typedef typename target_traits::box_iterator  target_box_iterator;
  //! Target tree body type
  typedef typename target_traits::body_type     target_body_type;
  //! Target tree body iterator
  typedef typename target_traits::body_iterator target_body_iterator;
};


#define FMMTL_IMPORT_TREEPAIR_TRAITS(TREETRAITS)                     \
  typedef typename TREETRAITS::source_tree_type     source_tree_type;   \
  typedef typename TREETRAITS::target_tree_type     target_tree_type;   \
  typedef typename TREETRAITS::source_box_type      source_box_type;    \
  typedef typename TREETRAITS::source_box_iterator  source_box_iterator; \
  typedef typename TREETRAITS::source_body_type     source_body_type;   \
  typedef typename TREETRAITS::source_body_iterator source_body_iterator; \
  typedef typename TREETRAITS::target_box_type      target_box_type;    \
  typedef typename TREETRAITS::target_box_iterator  target_box_iterator; \
  typedef typename TREETRAITS::target_body_type     target_body_type;   \
  typedef typename TREETRAITS::target_body_iterator target_body_iterator
