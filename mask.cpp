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

#include "mask.h"

FixedMask::FixedMask()
 : FixedImage<bool>(), _internal(0)
{

}


FixedMask::FixedMask(uint size_x, uint size_y)
 : FixedImage<bool>(size_x, size_y, false)
{
	init_internal(Point(size_x, size_y), Point(-1, -1), true, Point(-1, -1), Point(-2, -2), true);
}


FixedMask::FixedMask(uint size_x, uint size_y, bool default_value)
 : FixedImage<bool>(size_x, size_y, default_value)
{
	if (default_value) {
		init_internal(Point(0, 0), Point(size_x - 1, size_y - 1), true, Point(0, 0), Point(size_x - 1, size_y - 1), true);
	} else {
		init_internal(Point(size_x, size_y), Point(-1, -1), true, Point(-1, -1), Point(-2, -2), true);
	}
}


FixedMask::FixedMask(Shape size)
 : FixedImage<bool>(size, false)
{
	init_internal(Point(size.size_x, size.size_y), Point(-1, -1), true, Point(-1, -1), Point(-2, -2), true);
}


FixedMask::FixedMask(Shape size, bool default_value)
 : FixedImage<bool>(size.size_x, size.size_y, default_value)
{
	if (default_value) {
		init_internal(Point(0, 0), Point(size.size_x - 1, size.size_y - 1), true, Point(0, 0), Point(size.size_x - 1, size.size_y - 1), true);
	} else {
		init_internal(Point(size.size_x, size.size_y), Point(-1, -1), true, Point(-1, -1), Point(-2, -2), true);
	}
}


/**
 * Deep copy
 */
FixedMask::FixedMask(const Image<bool> &source)
 : FixedImage<bool>(source.clone())	// we have to invoke deep copying explicitly, because Image to FixedImage cast leads to the data sharing
{
	init_internal(Point::empty, Point::empty, false, Point::empty, Point::empty, false);
}


/**
 * Deep copy
 */
FixedMask::FixedMask(const FixedImage<bool> &source)
 : FixedImage<bool>(source.clone())	// we have to invoke deep copying explicitly, because FixedImage to FixedImage cast leads to the data sharing
{
	init_internal(Point::empty, Point::empty, false, Point::empty, Point::empty, false);
}


/**
 * Ref++
 */
FixedMask::FixedMask(const Mask &source)
 : FixedImage<bool>(source)	// Mask is casted to FixedImage implicitly, then FixedImage to FixedImage cast leads to the data sharing
{
	_internal = source._internal;
}


/**
 * Ref++
 */
FixedMask::FixedMask(const FixedMask &source)
 : FixedImage<bool>(source)	// FixedMask is casted to FixedImage implicitly, then FixedImage to FixedImage cast leads to the data sharing
{
	_internal = source._internal;
}


FixedMask::~FixedMask()
{
	// NOTE: _internal will be released in the destroy(), if needed
}


FixedMask& FixedMask::operator= (const Mask &other)
{
	// check for self-assignment
	if(this == &other) {
		return *this;
	}

	// finish all deals with the previous data
	release();

	// increase counter
	if (other._ref) {
		other._ref->counter++;
	}

	// assign new data
	this->_ref = other._ref;
	this->_size_x = other._size_x;
	this->_size_y = other._size_y;
	this->_number_of_channels = other._number_of_channels;
	this->_data = other._data;
	this->_internal = other._internal;

	return *this;
}


FixedMask& FixedMask::operator= (const FixedMask &other)
{
	// check for self-assignment
	if(this == &other) {
		return *this;
	}

	// finish all deals with the previous data
	release();

	// increase counter
	if (other._ref) {
		other._ref->counter++;
	}

	// assign new data
	this->_ref = other._ref;
	this->_size_x = other._size_x;
	this->_size_y = other._size_y;
	this->_number_of_channels = other._number_of_channels;
	this->_data = other._data;
	this->_internal = other._internal;

	return *this;
}


FixedMask::iterator FixedMask::begin() const
{
	return iterator(this, first());
}


FixedMask::iterator FixedMask::end() const
{
	return iterator(this, Point(_size_x, _size_y));
}


FixedMask::iterator FixedMask::rbegin() const
{
	return iterator(this, last(), true);
}


FixedMask::iterator FixedMask::rend() const
{
	return iterator(this, Point(-1, -1), true);
}


/**
 * Returns the element without range checking. Does not affect internal cache.
 */
bool FixedMask::get(uint x, uint y) const
{
	return _data[get_index(x, y, 0)];
}


/**
 * Returns the element without range checking. Does not affect internal cache.
 */
bool FixedMask::get(uint x, uint y, uint channel) const
{

	return _data[get_index(x, y, channel)];
}


/**
 * Returns the element without range checking. Does not affect internal cache.
 */
bool FixedMask::get(Point p) const
{
	return _data[get_index(p.x, p.y, 0)];
}


/**
 * Returns the element without range checking. Does not affect internal cache.
 */
bool FixedMask::get(Point p, uint channel) const
{
	return _data[get_index(p.x, p.y, channel)];
}


/**
 * Returns the element with range checking. Does not affect internal cache.
 * @return Value at a given coordinates or 'false', if out of range.
 */
bool FixedMask::test(uint x, uint y) const
{
	if (x >= _size_x || y >= _size_y) {
		return false;
	}

	return _data[get_index(x, y, 0)];
}


/**
 * Returns the element with range checking. Does not affect internal cache.
 * @return Value at a given coordinates or 'false', if out of range.
 */
bool FixedMask::test(uint x, uint y, uint channel) const
{
	if (x >= _size_x || y >= _size_y || channel >= _number_of_channels) {
		return false;
	}

	return _data[get_index(x, y, channel)];
}


/**
 * Returns the element with range checking. Does not affect internal cache.
 * @return Value at a given coordinates or 'false', if out of range.
 */
bool FixedMask::test(Point p) const
{
	if (p.x < 0 || (uint)p.x >= _size_x || p.y < 0 || (uint)p.y >= _size_y) {
		return false;
	}

	return _data[get_index(p.x, p.y, 0)];
}


/**
 * Returns the element with range checking. Does not affect internal cache.
 * @return Value at a given coordinates or 'false', if out of range.
 */
bool FixedMask::test(Point p, uint channel) const
{
	if (p.x < 0 || (uint)p.x >= _size_x || p.x < 0 || (uint)p.y >= _size_y || channel >= _number_of_channels) {
		return false;
	}

	return _data[get_index(p.x, p.y, channel)];
}


/**
 * Invokes deep copy.
 */
FixedMask FixedMask::clone() const
{
	FixedMask clone;
	clone._size_x = this->_size_x;
	clone._size_y = this->_size_y;
	clone._number_of_channels = this->_number_of_channels;

	if (this->_ref) {
		clone.init(this->_size_x, this->_size_y, this->_number_of_channels);
		memcpy(clone._data, this->_data,  this->_number_of_channels * this->_size_y * this->_size_x * sizeof(bool));
		clone.init_internal(_internal->first, _internal->last, _internal->is_first_last_valid,
							_internal->top_left, _internal->bottom_right, _internal->is_bounding_box_valid);
	}

	return clone;
}


/**
 * Invokes deep copy and then invert the cloned mask.
 */
FixedMask FixedMask::clone_invert() const
{
	FixedMask clone;
	clone._size_x = this->_size_x;
	clone._size_y = this->_size_y;
	clone._number_of_channels = this->_number_of_channels;

	if (this->_ref) {
		clone.init(this->_size_x, this->_size_y, this->_number_of_channels);
		for (uint i = 0; i < _size_x * _size_y * _number_of_channels; i++) {
			clone._data[i] = !this->_data[i];
		}

		clone.init_internal(Point::empty, Point::empty, false, Point::empty, Point::empty, false);
	}

	return clone;
}


/**
 * Returns masked points as a vector.
 */
vector<Point> FixedMask::get_masked_points() const
{
	if (!_internal->is_points_cache_valid) {
		actualize_points_cache();
	}

	return _internal->points_cache;
}


Point FixedMask::first() const
{
	if (!_internal->is_first_last_valid) {
		actualize_first_last();
	}

	return _internal->first;
}


Point FixedMask::last() const
{
	if (!_internal->is_first_last_valid) {
		actualize_first_last();
	}

	return _internal->last;
}


Point FixedMask::next(const Point &current) const
{
	int from_x = current.x + 1;
	for (uint y = current.y; y < _size_y; y++) {
		for (uint x = from_x; x < _size_x; x++) {
			if (_data[get_index(x, y, 0)]) {
				return Point(x, y);
			}
		}
		from_x = 0;
	}

	return Point(_size_x, _size_y);
}


Point FixedMask::prev(const Point &current) const
{
	int from_x = current.x - 1;
	for (int y = current.y; y >= 0; y--) {
		for (int x = from_x; x >= 0; x--) {
			if (_data[get_index(x, y, 0)]) {
				return Point(x, y);
			}
		}
		from_x = _size_x - 1;
	}

	return Point(-1, -1);
}


Point FixedMask::bounding_box_top_left() const
{
	if (!_internal->is_bounding_box_valid) {
		actualize_points_cache();
	}

	return _internal->top_left;
}


Point FixedMask::bounding_box_bottom_right() const
{
	if (!_internal->is_bounding_box_valid) {
		actualize_points_cache();
	}

	return _internal->bottom_right;
}


/* Protected */

void FixedMask::destroy() const
{
	delete _internal;
	FixedImage<bool>::destroy();
}

inline void FixedMask::init_internal(Point first, Point last, bool is_first_last_valid, Point top_left, Point bottom_rigth, bool is_bounding_box_valid)
{
	_internal = new __Internal();
	_internal->first = first;
	_internal->last = last;
	_internal->top_left = top_left;
	_internal->bottom_right = bottom_rigth;
	_internal->is_first_last_valid = is_first_last_valid;
	_internal->is_bounding_box_valid = is_bounding_box_valid;
	_internal->is_points_cache_valid = false;
}


inline void FixedMask::actualize_first_last() const
{
	_internal->first = Point(_size_x, _size_y);
	_internal->last = Point(-1, -1);

	bool is_last_found = false;
	for (int y = _size_y - 1; (!is_last_found) && (y >= 0); y--) {
		for (int x = _size_x - 1; (!is_last_found) && (x >= 0); x--) {
			if (_data[get_index(x, y, 0)]) {
				_internal->last = Point(x, y);
				is_last_found = true;
			}
		}
	}

	if (is_last_found) {
		bool _is_first_found = false;
		for (uint y = 0; (!_is_first_found) && (y < _size_y); y++) {
			for (uint x = 0; (!_is_first_found) && (x < _size_x); x++) {
				if (_data[get_index(x, y, 0)]) {
					_internal->first = Point(x, y);
					_is_first_found = true;
				}
			}
		}
	}

	_internal->is_first_last_valid = true;
}


/**
 * Populates cache by masked points and (as a side effect) updates bounding box
 */
inline void FixedMask::actualize_points_cache() const
{
	_internal->points_cache.clear();

	Point top_left(_size_x, _size_y);
	Point bottom_right(-1, -1);

	// fill points cache
	for(iterator it = begin(); it != end(); ++it) {
		_internal->points_cache.push_back(*it);

		top_left.x = min(top_left.x, it->x);
		top_left.y = min(top_left.y, it->y);
		bottom_right.x = max(bottom_right.x, it->x);
		bottom_right.y = max(bottom_right.y, it->y);
	}

	if (top_left <= bottom_right) {
		_internal->top_left = top_left;
		_internal->bottom_right = bottom_right;
	} else {
		_internal->top_left = Point(-1, -1);
		_internal->bottom_right = Point(-1, -1);
	}

	_internal->is_bounding_box_valid = true;
	_internal->is_points_cache_valid = true;
}


/* ==================== Mask ==================== */

Mask::Mask()
	: FixedMask()
{

}


Mask::Mask(uint size_x, uint size_y)
	: FixedMask(size_x, size_y, false)
{

}


Mask::Mask(uint size_x, uint size_y, bool default_value)
	: FixedMask(size_x, size_y, default_value)
{

}


Mask::Mask(Shape size)
	: FixedMask(size, false)
{

}


Mask::Mask(Shape size, bool default_value)
	: FixedMask(size.size_x, size.size_y, default_value)
{

}


/**
 * Deep copy
 */
Mask::Mask(const Image<bool> &source)
	: FixedMask(source)	// Image to FixedMask cast leads to deep copying
{

}


/**
 * Deep copy
 */
Mask::Mask(const FixedImage<bool> &source)
	: FixedMask(source)	// FixedImage to FixedMask cast leads to deep copying
{

}


/**
 * Ref++
 */
Mask::Mask(const Mask &source)
	: FixedMask(source)	// Mask to FixedMask cast leads to the data sharing
{

}


/**
 * Deep copy
 */
Mask::Mask(const FixedMask &source)
	: FixedMask(source.clone())	// we have to invoke deep copying explicitly, because FixedImage to FixedImage cast leads to the data sharing
{

}


Mask::~Mask()
{

}


Mask& Mask::operator= (const Mask &other)
{
	// check for self-assignment
	if(this == &other) {
		return *this;
	}

	// finish all deals with the previous data
	release();

	// increase counter
	if (other._ref) {
		other._ref->counter++;
	}

	// assign new data
	this->_ref = other._ref;
	this->_size_x = other._size_x;
	this->_size_y = other._size_y;
	this->_number_of_channels = other._number_of_channels;
	this->_data = other._data;
	this->_internal = other._internal;

	return *this;
}


Mask& Mask::operator= (const FixedMask &other)
{
	// check for self-assignment
	if(this == &other) {
		return *this;
	}

	// finish all deals with the previous data
	release();

	if (other._ref) {
		this->_size_x = other._size_x;
		this->_size_y = other._size_y;
		this->_number_of_channels = other._number_of_channels;
		init(other._size_x, other._size_y, other._number_of_channels);
		memcpy(this->_data, other._data,  other._number_of_channels * other._size_y * other._size_x * sizeof(bool));
		init_internal(other._internal->first, other._internal->last, other._internal->is_first_last_valid,
					  other._internal->top_left, other._internal->bottom_right, other._internal->is_bounding_box_valid);
	} else {
		this->_size_x = 0;
		this->_size_y = 0;
		this->_number_of_channels = 0;
		this->_ref = 0;
		this->_data = 0;
		this->_internal = 0;
	}

	return *this;
}


/**
 * Returns a reference to the element without range checking. Invalidates internal cache.
 */
bool& Mask::operator() (uint x, uint y)
{
	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	return _data[get_index(x, y, 0)];
}


/**
 * Returns a reference to the element without range checking. Invalidates internal cache.
 */
bool& Mask::operator() (uint x, uint y, uint channel)
{
	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	return _data[get_index(x, y, channel)];
}


/**
 * Returns a reference to the element without range checking. Invalidates internal cache.
 */
bool& Mask::operator() (Point p)
{
	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	return _data[get_index(p.x, p.y, 0)];
}


/**
 * Returns a reference to the element without range checking. Invalidates internal cache.
 */
bool& Mask::operator() (Point p, uint channel)
{
	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	return _data[get_index(p.x, p.y, channel)];
}


/**
 * Returns a reference to the element with range checking. Invalidates internal cache.
 * @note Throws std::out_of_range exception, if out of range.
 */
bool& Mask::at(uint x, uint y)
{
	if (x >= _size_x || y >= _size_y) {
		throw std::out_of_range("x or y coordinate is out of range");
	}

	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	return _data[get_index(x, y, 0)];
}


/**
 * Returns a reference to the element with range checking. Invalidates internal cache.
 * @note Throws std::out_of_range exception, if out of range.
 */
bool& Mask::at(uint x, uint y, uint channel)
{
	if (x >= _size_x || y >= _size_y || channel >= _number_of_channels) {
		throw std::out_of_range("channel, x or y coordinate is out of range");
	}

	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	return _data[get_index(x, y, channel)];
}


/**
 * Returns a reference to the element with range checking. Invalidates internal cache.
 * @note Throws std::out_of_range exception, if out of range.
 */
bool& Mask::at(Point p)
{
	if (p.x < 0 || (uint)p.x >= _size_x || p.y < 0 || (uint)p.y >= _size_y) {
		throw std::out_of_range("x or y coordinate is out of range");
	}

	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	return _data[get_index(p.x, p.y, 0)];
}


/**
 * Returns a reference to the element with range checking. Invalidates internal cache.
 * @note Throws std::out_of_range exception, if out of range.
 */
bool& Mask::at(Point p, uint channel)
{
	if (p.x < 0 || (uint)p.x >= _size_x || p.x < 0 || (uint)p.y >= _size_y || channel >= _number_of_channels) {
		throw std::out_of_range("channel, x or y coordinate is out of range");
	}

	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	return _data[get_index(p.x, p.y, channel)];
}


/**
 * Sets element to 'true' with range checking (does nothing, if out of range).
 * Invalidates internal cache.
 */
void Mask::mask(uint x, uint y)
{
	if (x >= _size_x || y >= _size_y) {
		return;
	}

	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	_data[get_index(x, y, 0)] = true;
}


/**
 * Sets element to 'true' with range checking (does nothing, if out of range).
 * Invalidates internal cache.
 */
void Mask::mask(uint x, uint y, uint channel)
{
	if (x >= _size_x || y >= _size_y || channel >= _number_of_channels) {
		return;
	}

	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	_data[get_index(x, y, channel)] = true;
}


/**
 * Sets element to 'true' with range checking (does nothing, if out of range).
 * Invalidates internal cache.
 */
void Mask::mask(Point p)
{
	if (p.x < 0 || (uint)p.x >= _size_x || p.y < 0 || (uint)p.y >= _size_y) {
		return;
	}

	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	_data[get_index(p.x, p.y, 0)] = true;
}


/**
 * Sets element to 'true' with range checking (does nothing, if out of range).
 * Invalidates internal cache.
 */
void Mask::mask(Point p, uint channel)
{
	if (p.x < 0 || (uint)p.x >= _size_x || p.x < 0 || (uint)p.y >= _size_y || channel >= _number_of_channels) {
		return;
	}

	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	_data[get_index(p.x, p.y, channel)] = true;
}


/**
 * Sets element to 'false' with range checking (does nothing, if out of range).
 * Invalidates internal cache.
 */
void Mask::unmask(uint x, uint y)
{
	if (x >= _size_x || y >= _size_y) {
		return;
	}

	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	_data[get_index(x, y, 0)] = false;
}


/**
 * Sets element to 'false' with range checking (does nothing, if out of range).
 * Invalidates internal cache.
 */
void Mask::unmask(uint x, uint y, uint channel)
{
	if (x >= _size_x || y >= _size_y || channel >= _number_of_channels) {
		return;
	}

	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	_data[get_index(x, y, channel)] = true;
}


/**
 * Sets element to 'false' with range checking (does nothing, if out of range).
 * Invalidates internal cache.
 */
void Mask::unmask(Point p)
{
	if (p.x < 0 || (uint)p.x >= _size_x || p.y < 0 || (uint)p.y >= _size_y) {
		return;
	}

	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	_data[get_index(p.x, p.y, 0)] = true;
}


/**
 * Sets element to 'false' with range checking (does nothing, if out of range).
 * Invalidates internal cache.
 */
void Mask::unmask(Point p, uint channel)
{
	if (p.x < 0 || (uint)p.x >= _size_x || p.x < 0 || (uint)p.y >= _size_y || channel >= _number_of_channels) {
		return;
	}

	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;

	_data[get_index(p.x, p.y, channel)] = true;
}


/**
 * Invokes deep copy.
 */
Mask Mask::clone() const
{
	Mask clone;
	clone._size_x = this->_size_x;
	clone._size_y = this->_size_y;
	clone._number_of_channels = this->_number_of_channels;

	if (this->_ref) {
		clone.init(this->_size_x, this->_size_y, this->_number_of_channels);
		memcpy(clone._data, this->_data,  this->_number_of_channels * this->_size_y * this->_size_x * sizeof(bool));

		clone.init_internal(_internal->first, _internal->last, _internal->is_first_last_valid,
							_internal->top_left, _internal->bottom_right, _internal->is_bounding_box_valid);
	}

	return clone;
}


/**
 * Invokes deep copy and then invert the cloned mask.
 * @return Inverted copy of current mask.
 */
Mask Mask::clone_invert() const
{
	Mask clone;
	clone._size_x = this->_size_x;
	clone._size_y = this->_size_y;
	clone._number_of_channels = this->_number_of_channels;

	if (this->_ref) {
		clone.init(this->_size_x, this->_size_y, this->_number_of_channels);
		for (uint i = 0; i < _size_x * _size_y * _number_of_channels; i++) {
			clone._data[i] = !this->_data[i];
		}

		clone.init_internal(Point::empty, Point::empty, false, Point::empty, Point::empty, false);
	}

	return clone;
}


/**
 * Invert current mask.
 */
void Mask::invert()
{
	for (uint i = 0; i < _size_x * _size_y * _number_of_channels; i++) {
		_data[i] = !_data[i];
	}

	_internal->is_first_last_valid = false;
	_internal->is_points_cache_valid = false;
}

/* Protected */

void Mask::destroy() const
{
	// no additional data
	FixedMask::destroy();
}
