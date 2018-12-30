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


#include "point.h"

Point Point::empty = Point(-1, -1);


Point::Point()
{
	x = 0;
	y = 0;
}


Point::Point(int x, int y)
{
	this->x = x;
	this->y = y;
}


bool Point::operator== (const Point &p) const
{
	return (this->x == p.x) && (this->y == p.y);
}


bool Point::operator!= (const Point &p) const
{
	return !((*this) == p);
}


Point& Point::operator= (const Point &p)
{
	if (this != &p) {
		this->x = p.x;
		this->y = p.y;
	}

	return *this;
}


Point& Point::operator+= (const Point &p)
{
	x += p.x;
	y += p.y;
	return *this;
}


Point& Point::operator-= (const Point &p)
{
	x -= p.x;
	y -= p.y;
	return *this;
}

const Point Point::operator+ (const Point &p) const
{
	Point result = *this;
	result += p;
	return result;
}


const Point Point::operator- (const Point &p) const
{
	Point result = *this;
	result -= p;
	return result;
}
