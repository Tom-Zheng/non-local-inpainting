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

#include "mask_iterator.h"

MaskIterator::MaskIterator()
{
	_mask = 0;
	_current = Point::empty;
	_is_reverse = false;
}


MaskIterator::MaskIterator(const FixedMask *mask, Point current, bool reverse)
{
	_mask = mask;
	_current = current;
	_is_reverse = reverse;
}


MaskIterator::MaskIterator(const MaskIterator& source)
{
	this->_mask = source._mask;
	this->_current = source._current;
	this->_is_reverse = source._is_reverse;
}


MaskIterator::~MaskIterator()
{

}


MaskIterator& MaskIterator::operator=(const MaskIterator& source)
{
	if (this != &source) {
		this->_mask = source._mask;
		this->_current = source._current;
		this->_is_reverse = source._is_reverse;
	}

	return *this;
}


bool MaskIterator::operator==(const MaskIterator& other) const
{
	return (this->_mask == other._mask) && (this->_current == other._current) && (this->_is_reverse == other._is_reverse);
}


bool MaskIterator::operator!=(const MaskIterator& other) const
{
	return (this->_mask != other._mask) || (this->_current != other._current) || (this->_is_reverse != other._is_reverse);
}


MaskIterator& MaskIterator::operator++()
{
	if (_mask) {
		_current = (_is_reverse) ?
					_mask->prev(_current) :
					_mask->next(_current);
	}

	return *this;
}


MaskIterator MaskIterator::operator++(int)
{
	MaskIterator aux(*this);

	if (_mask) {
		_current = (_is_reverse) ?
					_mask->prev(_current) :
					_mask->next(_current);
	}

	return aux;
}


MaskIterator& MaskIterator::operator--()
{
	if (_mask) {
		_current = (_is_reverse) ?
					_mask->next(_current) :
					_mask->prev(_current);
	}

	return *this;
}


MaskIterator MaskIterator::operator--(int)
{
	MaskIterator aux(*this);

	if (_mask) {
		_current = (_is_reverse) ?
					_mask->next(_current) :
					_mask->prev(_current);
	}

	return aux;
}


const MaskIterator::const_reference MaskIterator::operator*() const
{
	return _current;
}


const MaskIterator::const_pointer MaskIterator::operator->() const
{
	return &_current;
}
