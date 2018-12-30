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

#include "shape.h"

Shape Shape::empty = Shape();

Shape::Shape()
{
	size_x = 0;
	size_y = 0;
}


Shape::Shape(unsigned int size_x, unsigned int size_y)
{
	this->size_x = size_x;
	this->size_y = size_y;
}


Shape::Shape(const Shape &source)
{
	size_x = source.size_x;
	size_y = source.size_y;
}


bool Shape::operator== (const Shape &other) const
{
	return (this->size_x == other.size_x) && (this->size_y == other.size_y);
}


bool Shape::operator!= (const Shape &other) const
{
	return !((*this) == other);
}


bool Shape::is_empty() const
{
	return size_x <= 0 || size_y <= 0;
}


/**
 * True if point lies inside the shape (relative to the top left point).
 */
bool Shape::contains(const Point &p) const
{
	return p.x >= 0 && p.y >= 0 && (unsigned)p.x < size_x && (unsigned)p.y < size_y;
}


/**
 * True if point lies inside the shape (relative to the top left point).
 */
bool Shape::contains(int x, int y) const
{
	return x >= 0 && y >= 0 && (unsigned)x < size_x && (unsigned)y < size_y;
}


/**
 * True if point lies inside the shape (relative to the center point).
 */
bool Shape::abs_contains(const Point &p) const
{
	return std::abs((float) p.x) * 2 < size_x && std::abs((float) p.y) * 2 < size_y;
}



