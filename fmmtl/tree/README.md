Spacial Tree Datatypes
=====
Trees are immutable data structures representing hierarchical paritions of a set of points. The trees in this library do not actually store the points or any corresponding data, but encapsulate only a permutation of input range of "bodies" it was constructed on and the hierarchically disjoint subsets of the permuted range representing "boxes".

Formally, a tree on a dataset of "bodies" ![equation](http://latex.codecogs.com/gif.latex?%5Cmathcal%7BP%7D%20%3D%20%5C%7Bp_0%2C%20p_1%2C%20%5Cldots%2C%20p_%7Bn-1%7D%5C%7D) is a directed, connected, acyclic, rooted simple graph with the following proerties:
* Each node (or "box" in the implementation), holds a zero or more bodies, is connected to exactly one parent node, and has zero or more children.
* The node with itself as the parent node is the root of the tree.
* Each body in ![equation](http://latex.codecogs.com/gif.latex?%5Cmathcal%7BP%7D) is contained in the root.
* The child nodes contain mutually disjoint subsets of the set of bodies contained by the parent.

The tree concept is
```c++
concept Tree {
  typedef size_type;
  typedef box_type;
  typedef body_type;
  typedef point_type;

  typedef box_iterator;
  typedef body_iterator;
  typedef permute_iterator;

  // Constructors...

  point_type              center() const;
  BoundingBox<point_type> extents() const;

  size_type size() const;
  size_type bodies() const;
  size_type boxes() const;
  size_type levels() const;

  box_type     root() const;
  box_iterator box_begin() const;
  box_iterator box_end() const;
  box_iterator box_begin(size_type level) const;
  box_iterator box_end(size_type level) const;

  body_iterator body_begin() const;
  body_iterator body_end() const;

  permute_iterator permute_begin() const;
  permute_iterator permute_end() const;
};
```
where the `box_type` has the concept
```c++
concept box_type {
  point_type              center() const;
  point_type::value_type  radius_sq() const;
  BoundingBox<point_type> extents() const;

  box_type      parent() const;
  box_iterator  child_begin() const;
  box_iterator  child_end() const;
  body_iterator body_begin() const;
  body_iterator body_end() const;

  bool operator==(const box_type&) const;
  bool operator <(const box_type&) const;
}
```
and the `body_type` has the concept
```c++
concept body_type {
  size_type index() const;
  size_type number() const;
}
```
