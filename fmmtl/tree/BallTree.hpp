#pragma once

#include <vector>
#include <algorithm>
#include <iterator>
#include <tuple>
#include <iostream>

#include <boost/range.hpp>
#include <boost/iterator/iterator_adaptor.hpp>
#include <boost/iterator/counting_iterator.hpp>

#include "fmmtl/util/Logger.hpp"
#include "fmmtl/numeric/Vec.hpp"
#include "fmmtl/tree/BoundingSphere.hpp"
#include "fmmtl/numeric/bits.hpp"
#include "fmmtl/numeric/norm.hpp"

namespace fmmtl {
using boost::has_range_iterator;
using boost::iterator_adaptor;
using boost::counting_iterator;

template <unsigned DIM>
class BallTree {
	struct Box;
	struct Body;
	struct BoxData;
	struct BoxIterator;
	struct BodyIterator;

	typedef BallTree<DIM> tree_type;

public:
	typedef unsigned size_type;
	typedef Vec<DIM, double> point_type;

	typedef Box box_type;
	typedef Body body_type;
	typedef BoxIterator box_iterator;
	typedef BodyIterator body_iterator;

private:
	//Representation
	std::vector<size_type> permute_;
	std::vector<BoxData> box_data_;

	struct BoxData {
		typedef BoundingSphere<point_type> bounding_sphere_type;
		//first body index in the box
		size_type body_begin_;
		//one-past-last body index in the box
		size_type body_end_;
		//Bounding sphere
		bounding_sphere_type bounding_sphere_;

		//memberwise constructor
		BoxData(size_type bb, size_type be, const bounding_sphere_type& bs)
			: body_begin_(bb), body_end_(be), bounding_sphere_(bs) {}

	};

	struct Body {
	public:
		//Construct an invalid Body
		Body() {}
		//original order of the body
		size_type number() const {
			return tree_->permute_[idx_];
		}
		//current order
		size_type index() const {
			return idx_;
		}

	private:
		size_type idx_;
		tree_type* tree_;
		Body(size_type idx, const tree_type* tree)
			: idx_(idx), tree_(const_cast<tree_type*>(tree)) {
			FMMTL_ASSERT(idx_<tree_->size());
		}
		friend class BallTree;
	};

	struct Box {
	public:
		typedef typename tree_type::box_iterator box_iterator;
		typedef typename tree_type::body_iterator body_iterator;

		//Construct an invalid Box
		Box() {}

		size_type index() const {
			return idx_;
		}
		size_type level() const {
			return std::log2(idx_+1);
		}
		double volume() const {
			//XXX: need this function?
			return 0;
		}

		//return the radius of this bounding sphere
		double radius() const {
			return data().bounding_sphere_.radius();
		}
		//return the center coordinate of this bounding sphere
		point_type center() const {
			return data().bounding_sphere_.center();
		}
		Box parent() const {
			FMMTL_ASSERT(!(*this == tree_->root()));
			return Box ((idx_ - 1)/2, tree_);
		}
		bool is_leaf() const {
			return (2*idx_+1 >= tree_->boxes());
		}

		box_iterator child_begin() const {
			FMMTL_ASSERT(!is_leaf());
			return box_iterator(2*idx_+1, tree_);
		}
		box_iterator child_end() const {
			FMMTL_ASSERT(!is_leaf());
			return box_iterator(2*idx_+3, tree_);
		}
		body_iterator body_begin() const {
			return body_iterator(data().body_begin_, tree_);
		}
		body_iterator body_end() const {
			return body_iterator(data().body_end_, tree_);
		}

		constexpr size_type num_children() const {
			return 2;	//XXX: box in a full binary tree always has 2 children
		}
		size_type num_bodies() const {
			return std::distance(body_begin(), body_end());
		}

		bool operator==(const Box& rhs) const {
			FMMTL_ASSERT(tree_==rhs.tree_);
			return idx_ == rhs.idx_;
		}
		bool operator<(const Box& rhs) const {
			FMMTL_ASSERT(tree_==rhs.tree_);
			return idx_ < rhs.idx_;
		}

		inline friend std::ostream& operator<<(std::ostream& s, const box_type& b) {
			size_type num_bodies = b.num_bodies();
			size_type first_body = b.body_begin()->index();
			size_type last_body = first_body + num_bodies - 1;
			size_type parent_idx = b.index()==0 ? 0 : b.parent().index();

			return s << "Box " << b.index()
					 << " (L" << b.level() << ", P" << parent_idx
					 << ", " << num_bodies << (num_bodies == 1 ? " body" : " bodies")
					 << " " << first_body << "-" << last_body
					 << "): Center: " << b.center() << ";\t Radius: " << b.radius();
		}

	private:
		size_type idx_;
		tree_type* tree_;
		Box(size_type idx, const tree_type* tree)
			: idx_(idx), tree_(const_cast<tree_type*>(tree)) {
			FMMTL_ASSERT(idx_<tree_->boxes());
		}
		inline BoxData& data() const {
			return tree_->box_data_[idx_];
		}
		friend class BallTree;
	};	//end struct Box

	struct BoxIterator
		: public iterator_adaptor<	BoxIterator,						//Derived class
									counting_iterator<size_type>,		//BaseType
									Box,								//Value
									std::random_access_iterator_tag,	//iterator_category
									Box>								//Reference
	{
	public:
		inline BoxIterator() {}
		inline size_type index() const {
			return *(this->base_reference());
		}

	private: 
		const tree_type* tree_;
		friend class BallTree;
		inline BoxIterator(size_type idx, const tree_type* tree)
			: BoxIterator::iterator_adaptor(counting_iterator<size_type>(idx)), tree_(tree) {}
		inline BoxIterator(const Box& b)
			: BoxIterator::iterator_adaptor(counting_iterator<size_type>(b.idx_)), tree_(b.tree_) {}

		friend class boost::iterator_core_access;
		inline Box dereference() const {
			return Box(index(), tree_);
		}
	};

	struct BodyIterator 
		: public iterator_adaptor<	BodyIterator,						//Derived class
									counting_iterator<size_type>,		//BaseType
									Body,								//Value
									std::random_access_iterator_tag,	//iterator_category
									Body>								//Reference
	{
	public:
		inline BodyIterator() {}
		inline size_type index() const {
			return *(this->base_reference());
		}

	private: 
		const tree_type* tree_;
		friend class BallTree;
		inline BodyIterator(size_type idx, const tree_type* tree)
			: BodyIterator::iterator_adaptor(counting_iterator<size_type>(idx)), tree_(tree) {}
		inline BodyIterator(const Box& b)
			: BodyIterator::iterator_adaptor(counting_iterator<size_type>(b.idx_)), tree_(b.tree_) {}

		friend class boost::iterator_core_access;
		inline Body dereference() const {
			return Body(index(), tree_);
		}
	};

public:
	
	//Tree constructor with Range
	template<typename Range>
	BallTree(const Range& rng, size_type n_crit = 256, typename std::enable_if<has_range_iterator<Range>::value>::type* = 0)
		: BallTree(rng.begin(), rng.end(), n_crit) {}
	
	template<typename PointIter>
	BallTree(PointIter first, PointIter last, size_type n_crit = 256) {
		insert(first, last, n_crit);
	}

	//return root Bounding Sphere
	BoundingSphere<point_type> bounding_sphere() const {
		return box_data_[0].bounding_sphere_;
	}
	point_type center() const {
		return root().center();
	}

	inline size_type size() const {
		return permute_.size();
	}
	inline size_type bodies() const {
		return size();
	}
	inline size_type boxes() const {
		return box_data_.size();
	}
	//return the number of boxes on level L
	inline size_type boxes(size_type L) const {
		return (1 << L);	// 1<<L == pow(2, L), assuming it's a full binary tree
	}
	inline size_type levels() const {
		return std::log2(boxes());
	}
	//maximum possible level of any box in this tree
	inline static size_type max_level() {
		return size_type(-1);	//XXX: meaning?
	}
	inline bool contains(const box_type& box) const {
		return this==box.tree_;
	}
	inline bool contains(const body_type& body) const {
		return this==body.tree_;
	}

	//return the root box of this tree
	box_type root() const {
		return Box(0, this);
	}
	box_type box(size_type idx) const {
		FMMTL_ASSERT(idx < box_data_.size());
		return Box(idx, this);
	}
	body_type body(size_type idx) const {
		FMMTL_ASSERT(idx < size());
		return Body(idx, this);
	}

	body_iterator body_begin() const {
		return body_iterator(0, this);
	}
	body_iterator body_end() const {
		return body_iterator(bodies(), this);
	}

	box_iterator box_begin() const {
		return box_iterator(0, this);
	}
	box_iterator box_end() const {
		return box_iterator(boxes(), this);
	}
	//return an iterator to the first box on level L 
	box_iterator box_begin(size_type L) const {}
	//return an iterator to one past the last box on level L 
	box_iterator box_end(size_type L) const {}

	inline friend std::ostream& operator<<(std::ostream& s, const tree_type& t)
	{
		struct {
			inline std::ostream& print(std::ostream& s, const box_type& box) {
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
	template <typename PointIter>
	void insert(PointIter p_first, PointIter p_last, size_type NCRIT) {
		FMMTL_LOG("BallTree insert");
		//possibly convert N-dimensional point to polar coordinates and store them
		assert (p_first != p_last);

		// Create a point-idx pair vector
		typedef typename std::iterator_traits<PointIter>::value_type point_i_type;
		typedef std::pair<point_i_type, unsigned> point_t;

		std::vector<point_t> point;
		std::vector<unsigned> pivot_arry;
		
		if (std::is_same<typename std::iterator_traits<PointIter>::iterator_category,
			std::random_access_iterator_tag>::value) 
				point.reserve(std::distance(p_first, p_last));

		//construct root bounding sphere, with point iterators
		unsigned idx = 0;
		for (auto pi = p_first; pi!=p_last; ++pi, ++idx)
			point.emplace_back(*pi, idx);

		permute_.reserve(point.size());
		pivot_arry.reserve(point.size());

		unsigned leaves = ceil_pow_2((point.size()+NCRIT-1)/NCRIT);
		unsigned levels = std::log2(leaves);

		auto bs_params = calc_bs_params(point.begin(), point.end());	//0: center, 1: radius_sq, 2: pivot
		auto curr_bs = BoundingSphere<point_i_type>(std::get<0>(bs_params), std::get<1>(bs_params));
		pivot_arry.emplace_back(std::get<2>(bs_params));

		box_data_.reserve(2*leaves-1);
		box_data_.emplace_back(0, size_type(point.size()), curr_bs);

		for (unsigned k=0; k<(1<<levels)-1; ++k) {
			auto pivot = pivot_arry[k];

			//get the dimension with the largest spread (pivot)
			auto myless_comp = [pivot] (const point_t& a, const point_t& b) {
				return a.first[pivot] < b.first[pivot];
			};

			//partition the points in the dimension @a pivot, with the median
			auto p_begin = point.begin() + box_data_[k].body_begin_;
			auto p_end = point.begin() + box_data_[k].body_end_;
			auto p_mid = p_begin + (p_end - p_begin + 1)/2;
			std::nth_element(p_begin, p_mid, p_end, myless_comp);

			//build 2 child boxes
			unsigned mid = p_mid - point.begin();
			bs_params = calc_bs_params(p_begin, p_mid);	//0: center, 1: radius_sq, 2: pivot
			curr_bs = BoundingSphere<point_i_type>(std::get<0>(bs_params), std::get<1>(bs_params));
			pivot_arry.emplace_back(std::get<2>(bs_params));
			box_data_.emplace_back(box_data_[k].body_begin_, mid, curr_bs);

			bs_params = calc_bs_params(p_mid, p_end);	//0: center, 1: radius_sq, 2: pivot
			curr_bs = BoundingSphere<point_i_type>(std::get<0>(bs_params), std::get<1>(bs_params));
			pivot_arry.emplace_back(std::get<2>(bs_params));
			box_data_.emplace_back(mid, box_data_[k].body_end_, curr_bs);
		}

		std::cout << box_data_.size() << "\t" << 2*leaves-1 << std::endl;
		assert (box_data_.size() <= 2*leaves-1);

		for (auto&& p : point)
			permute_.push_back(p.second);
	}

	static unsigned max_dim(const point_type& p) {
		return std::max_element(p.begin(), p.end()-p.begin());
	}

	/** calculate bounding sphere parameters : center, radius_sq and pivot for a given range of points [first, last)
	 * @tparam PointIndexIter takes the form of std::vector<point_type, unsigned>::iterator
	 */
	template <typename PointIndexIter>
	std::tuple<point_type, double, unsigned> calc_bs_params(PointIndexIter first, PointIndexIter last) {
		unsigned n = std::distance(first, last);
		unsigned dim = (*first).first.size();
		//calculate maximum spread in each dimension
		auto min_pos((*first).first);	//stores minimum point position in each dimension
		auto max_pos(min_pos);	//stores maxmium point position in each dimension
		auto center = point_type(0,0,0);
		double radius_sq = 0.;

		//Update center for the BoundingSphere
		for (auto new_first=first; new_first!=last; ++new_first) {
			for (unsigned i=0; i<dim; ++i) {
				center[i] += (*new_first).first[i];
				if ((*new_first).first[i] < min_pos[i])
					min_pos[i] = (*new_first).first[i];
				if ((*new_first).first[i] > max_pos[i])
					max_pos[i] = (*new_first).first[i];
			}
		}
		for (unsigned i=0; i<dim; ++i)
			center[i] /= n;

		//calculate pivot
		auto temp_pos = max_pos-min_pos;
		unsigned pivot = std::max_element(temp_pos.begin(), temp_pos.end()) - temp_pos.begin();

		//Update radius_sq
		for (; first!=last; ++first) {
			double dist_sq = norm_2_sq(center - (*first).first);	
			if (dist_sq > radius_sq)
				radius_sq = dist_sq;
		}
		return std::make_tuple(center, radius_sq, pivot);
	}

	BallTree (const BallTree&) {};
	void operator=(const BallTree&) {};
};	//end BallTree

}	//end namespace fmmtl
