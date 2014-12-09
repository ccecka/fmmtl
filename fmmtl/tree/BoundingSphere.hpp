#include <algorithm>
#include <cmath>
#include <iostream>
#include "fmmtl/numeric/norm.hpp"
#include "fmmtl/meta/dimension.hpp"

namespace fmmtl {

template <typename POINT>
class BoundingSphere {
public:
	typedef POINT point_type;

	//invalid constructor
	BoundingSphere() {};

	//Constuructor
	template <typename PointIter>
	BoundingSphere(PointIter first, PointIter last)
		: center_(point_type()), radius_(0.), pivot_(0) {
		assert(first != last);
		insert(first, last);
	}

	const point_type& center() const {
		return center_;
	}

	const double& radius() const {
		return radius_;
	}

	const unsigned& pivot() const {
		return pivot_;
	}

	bool constains(const point_type& p) const {
		if (norm_2(p - center()) < radius()) return true;
		else return false;
	}

	//XXX overloaded insert function with Range
	template <typename PointIter>
	BoundingSphere& insert(PointIter first, PointIter last) {
		unsigned n = std::distance(first, last);
		unsigned dim = (*first).size();	//XXX static or constexpr member
		//calculate maximum spread in each dimension
		point_type min_pos(*first);	//stores minimum point position in each dimension
		point_type max_pos(*first);	//stores maxmium point position in each dimension

		//Update center_ for the BoundingSphere
		for (auto new_first=first; new_first!=last; ++new_first) {
			for (unsigned i=0; i<dim; ++i) {
				center_[i] += (*new_first)[i];
				if ((*new_first)[i] < min_pos[i])
					min_pos[i] = (*new_first)[i];
				if ((*new_first)[i] > max_pos[i])
					max_pos[i] = (*new_first)[i];
			}
		}
		for (unsigned i=0; i<dim; ++i)
			center_[i] /= n;

		//calculate pivot
		auto temp_pos = max_pos-min_pos;
		pivot_ = std::max_element(temp_pos.begin(), temp_pos.end()) - temp_pos.begin();

		//Update radius_
		for (; first!=last; ++first) {
			double dist = norm_2(center_ - *first);	
			if (dist > radius_)
				radius_ = dist;
		}
		return *this;
	}

private:
	point_type center_;
	double radius_;
	unsigned pivot_;

};//end class RoundBoundingSphere


template <typename P>
inline std::ostream& operator<<(std::ostream& s, const BoundingSphere<P>& b) {
	const unsigned dim = b.center().size();
	s << "[ Center: " << b.center()[0];
	for (unsigned i = 1; i != dim; ++i)
		s << ", " << b.center()[i];
	s << " ;\tRadius: " << b.radius();
	return s << ';';
}

}//end namespace fmmtl
