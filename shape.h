/**
 * Copyright (C) 2015, Vadim Fedorov <vadim.fedorov@upf.edu>
 * Copyright (C) 2015, Gabriele Facciolo <facciolo@ens-cachan.fr>
 * Copyright (C) 2015, Pablo Arias <pablo.arias@cmla.ens-cachan.fr>
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the simplified BSD
 * License. You should have received a copy of this license along
 * this program. If not, see
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#ifndef SHAPE_H_
#define SHAPE_H_

#include <cmath>
#include "point.h"

/**
 * Container for 2d shape (size_x is width, size_y is height).
 */
struct Shape
{
	unsigned int size_x, size_y;

	Shape();
	Shape(unsigned int size_x, unsigned int size_y);
	Shape(const Shape &source);

	bool operator== (const Shape &other) const;
	bool operator!= (const Shape &other) const;

	friend inline bool operator< (const Shape& lhs, const Shape& rhs);
	friend inline bool operator> (const Shape& lhs, const Shape& rhs);
	friend inline bool operator<=(const Shape& lhs, const Shape& rhs);
	friend inline bool operator>=(const Shape& lhs, const Shape& rhs);

	bool is_empty() const;
	bool contains(const Point &p) const;
	bool contains(int x, int y) const;
	bool abs_contains(const Point &p) const;

	static Shape empty;
};

// NOTE: definitions are in header in order to overload two argument versions.
inline bool operator< (const Shape& lhs, const Shape& rhs)
{
	return lhs.size_y < rhs.size_y ||
		(lhs.size_y == rhs.size_y && lhs.size_x < rhs.size_x);
}
inline bool operator> (const Shape& lhs, const Shape& rhs) { return operator< (rhs,lhs); }
inline bool operator<= (const Shape& lhs, const Shape& rhs) { return !operator> (lhs,rhs); }
inline bool operator>= (const Shape& lhs, const Shape& rhs) { return !operator< (lhs,rhs); }


#endif /* SHAPE_H_ */
