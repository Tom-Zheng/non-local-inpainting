/**
 * Copyright (C) 2015, Vadim Fedorov <vadim.fedorov@upf.edu>
 * Copyright (C) 2015, Gabriele Facciolo <facciolo@ens-cachan.fr>
 * Copyright (C) 2015, Pablo Arias <pablo.arias@upf.edu>
 *
 * This program is free software: you can use, modify and/or
 * redistribute it under the terms of the simplified BSD
 * License. You should have received a copy of this license along
 * this program. If not, see
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

#include "image.h"

/* ==================== FixedImage ==================== */


template <class T>
FixedImage<T>::FixedImage()
 : _size_x(0), _size_y(0), _number_of_channels(0), _data(0), _ref(0)
{
}


template <class T>
FixedImage<T>::FixedImage(uint size_x, uint size_y)
 : _size_x(size_x), _size_y(size_y), _number_of_channels(1)
{
	init(size_x, size_y, 1);
}


template <class T>
FixedImage<T>::FixedImage(uint size_x, uint size_y, uint number_of_channels)
 : _size_x(size_x), _size_y(size_y), _number_of_channels(number_of_channels)
{
	init(size_x, size_y, number_of_channels);
}


template <class T>
FixedImage<T>::FixedImage(uint size_x, uint size_y, T default_value)
 : _size_x(size_x), _size_y(size_y), _number_of_channels(1)
{
	init(size_x, size_y, 1);
	fill_internal(default_value);
}


template <class T>
FixedImage<T>::FixedImage(uint size_x, uint size_y, uint number_of_channels, T default_value)
 : _size_x(size_x), _size_y(size_y), _number_of_channels(number_of_channels)
{
	init(size_x, size_y, number_of_channels);
	fill_internal(default_value);
}


template <class T>
FixedImage<T>::FixedImage(Shape size)
 : _size_x(size.size_x), _size_y(size.size_y), _number_of_channels(1)
{
	init(size.size_x, size.size_y, 1);
}


template <class T>
FixedImage<T>::FixedImage(Shape size, uint number_of_channels)
 : _size_x(size.size_x), _size_y(size.size_y), _number_of_channels(number_of_channels)
{
	init(size.size_x, size.size_y, number_of_channels);
}


template <class T>
FixedImage<T>::FixedImage(Shape size, T default_value)
 : _size_x(size.size_x), _size_y(size.size_y), _number_of_channels(1)
{
	init(size.size_x, size.size_y, 1);
	fill_internal(default_value);
}


template <class T>
FixedImage<T>::FixedImage(Shape size, uint number_of_channels, T default_value)
 : _size_x(size.size_x), _size_y(size.size_y), _number_of_channels(number_of_channels)
{
	init(size.size_x, size.size_y, number_of_channels);
	fill_internal(default_value);
}


template <class T>
FixedImage<T>::FixedImage(const FixedImage<T> &source)
 : _size_x(source._size_x), _size_y(source._size_y), _number_of_channels(source._number_of_channels), _data(source._data), _ref(source._ref)
{
	if (_ref) {
		#pragma omp atomic
		_ref->counter++;
	}
}


template <class T>
FixedImage<T>::FixedImage(const Image<T> &source)
 : _size_x(source._size_x), _size_y(source._size_y), _number_of_channels(source._number_of_channels), _data(source._data), _ref(source._ref)
{
	if (_ref) {
		#pragma omp atomic
		_ref->counter++;
	}
}


template <class T>
FixedImage<T>::~FixedImage()
{
	release();
}


template <class T>
FixedImage<T>& FixedImage<T>::operator= (const FixedImage<T> &other)
{
	// check for self-assignment
	if(this == &other) {
		return *this;
	}

	// finish all deals with the previous data
	FixedImage<T>::release();

	// increase counter
	if (other._ref) {
		#pragma omp atomic
		other._ref->counter++;
	}

	// assign new data
	this->_ref = other._ref;
	this->_size_x = other._size_x;
	this->_size_y = other._size_y;
	this->_number_of_channels = other._number_of_channels;
	this->_data = other._data;

	return *this;
}


template <class T>
FixedImage<T>& FixedImage<T>::operator= (const Image<T> &other)
{
	// finish all deals with the previous data
	FixedImage<T>::release();

	// increase counter
	if (other._ref) {
		#pragma omp atomic
		other._ref->counter++;
	}

	// assign new data
	this->_ref = other._ref;
	this->_size_x = other._size_x;
	this->_size_y = other._size_y;
	this->_number_of_channels = other._number_of_channels;
	this->_data = other._data;

	return *this;
}


template <class T>
bool FixedImage<T>::operator== (const FixedImage<T> &other) const
{
	return this->_ref == other._ref;
}


template <class T>
bool FixedImage<T>::operator== (const Image<T> &other) const
{
	return this->_ref == other._ref;
}


/**
 * Checks if current image is empty.
 * @return Is the image empty?
 */
template <class T>
bool FixedImage<T>::is_empty() const
{
	return _ref == 0;
}


/**
 * Checks if current image is not empty.
 * @return Is the image not empty?
 */
template <class T>
bool FixedImage<T>::is_not_empty() const
{
	return _ref != 0;
}


template <class T>
uint FixedImage<T>::get_size_x() const
{
	return _size_x;
}


template <class T>
uint FixedImage<T>::get_size_y() const
{
	return _size_y;
}


template <class T>
Shape FixedImage<T>::get_size() const
{
	return Shape(_size_x, _size_y);
}


template <class T>
uint FixedImage<T>::get_number_of_channels() const
{
	return _number_of_channels;
}


/**
 * Returns read-only value without range checking.
 */
template <class T>
const T& FixedImage<T>::operator() (uint x, uint y) const
{
	return _data[get_index(x, y, 0)];
}


/**
 * Returns read-only value without range checking.
 */
template <class T>
const T& FixedImage<T>::operator() (uint x, uint y, uint channel) const
{
	return _data[get_index(x, y, channel)];
}


/**
 * Returns read-only value without range checking.
 */
template <class T>
const T& FixedImage<T>::operator() (Point p) const
{
	return _data[get_index(p.x, p.y, 0)];
}


/**
 * Returns read-only value without range checking.
 */
template <class T>
const T& FixedImage<T>::operator() (Point p, uint channel) const
{
	return _data[get_index(p.x, p.y, channel)];
}


/**
 * Returns read-only value with range checking.
 * @note Throws std::out_of_range exception, if out of range.
 */
template <class T>
const T& FixedImage<T>::at(uint x, uint y) const
{
	if (x >= _size_x || y >= _size_y) {
		throw std::out_of_range("x or y coordinate is out of range");
	}

	return _data[get_index(x, y, 0)];
}


/**
 * Returns read-only value with range checking.
 * @note Throws std::out_of_range exception, if out of range.
 */
template <class T>
const T& FixedImage<T>::at(uint x, uint y, uint channel) const
{
	if (x >= _size_x || y >= _size_y || channel >= _number_of_channels) {
		throw std::out_of_range("channel, x or y coordinate is out of range");
	}

	return _data[get_index(x, y, channel)];
}


/**
 * Returns read-only value with range checking.
 * @note Throws std::out_of_range exception, if out of range.
 */
template <class T>
const T& FixedImage<T>::at(Point p) const
{
	if (p.x >= _size_x || p.y >= _size_y) {
		throw std::out_of_range("x or y coordinate is out of range");
	}

	return _data[get_index(p.x, p.y, 0)];
}


/**
 * Returns read-only value with range checking.
 * @note Throws std::out_of_range exception, if out of range.
 */
template <class T>
const T& FixedImage<T>::at(Point p, uint channel) const
{
	if (p.x >= _size_x || p.y >= _size_y || channel >= _number_of_channels) {
		throw std::out_of_range("channel, x or y coordinate is out of range");
	}

	return _data[get_index(p.x, p.y, channel)];
}


/**
 * Checks if all indexes are in range and modifies the value parameter.
 * @param value [out] Value of an element at a given coordinates, if in range.
 * @return Are given coordinates in range?
 */
template <class T>
bool FixedImage<T>::try_get_value(uint x, uint y, T& value) const
{
	if (x >= _size_x || y >= _size_y || !_data) {
		return false;
	}

	value = _data[get_index(x, y, 0)];

	return true;
}


/**
 * Checks if all indexes are in range and modifies the value parameter.
 * @param value [out] Value of an element at a given coordinates, if in range.
 * @return Are given coordinates in range?
 */
template <class T>
bool FixedImage<T>::try_get_value(uint x, uint y, uint channel, T& value) const
{
	if (x >= _size_x || y >= _size_y || channel >= _number_of_channels || !_data) {
		return false;
	}

	value = _data[get_index(x, y, channel)];

	return true;
}


/**
 * Checks if all indexes are in range and modifies the value parameter.
 * @param value [out] Value of an element at a given coordinates, if in range.
 * @return Are given coordinates in range?
 */
template <class T>
bool FixedImage<T>::try_get_value(Point p, T& value) const
{
	if (p.x >= _size_x || p.y >= _size_y || !_data) {
		return false;
	}

	value = _data[get_index(p.x, p.y, 0)];

	return true;
}


/**
 * Checks if all indexes are in range and modifies the value parameter.
 * @param value [out] Value of an element at a given coordinates, if in range.
 * @return Are given coordinates in range?
 */
template <class T>
bool FixedImage<T>::try_get_value(Point p, uint channel, T& value) const
{
	if (p.x >= _size_x || p.y >= _size_y || channel >= _number_of_channels || !_data) {
		return false;
	}

	value = _data[get_index(p.x, p.y, channel)];

	return true;
}


/**
 * Returns pointer to internal data.
 */
template <class T>
const T* FixedImage<T>::raw() const
{
	return this->_data;
}


/**
 * Invokes deep copy.
 */
template <class T>
FixedImage<T> FixedImage<T>::clone() const
{
	FixedImage<T> clone;
	clone._size_x = this->_size_x;
	clone._size_y = this->_size_y;
	clone._number_of_channels = this->_number_of_channels;

	if (this->_ref) {
		clone.init(this->_size_x, this->_size_y, this->_number_of_channels);
		memcpy(clone._data, this->_data,  this->_number_of_channels * this->_size_y * this->_size_x * sizeof(T));
	}

	return clone;
}

/* Protected */

template <class T>
inline void FixedImage<T>::init(uint size_x, uint size_y, uint number_of_channels)
{
	_data = new T[size_x * size_y * number_of_channels];
	_ref = new __Ref();
	_ref->counter = 1;
}


template <class T>
void FixedImage<T>::fill_internal(const T &value)
{
	std::fill_n(_data, _number_of_channels * _size_y * _size_x, value);
}


template <class T>
void FixedImage<T>::release() const
{
	if (_ref) {
		#pragma omp critical(IMAGE_REF_COUNTING)
		{
			#pragma omp atomic
			_ref->counter--;
			if (_ref->counter == 0) {
				destroy();
			}
		}
	}
}


template <class T>
void FixedImage<T>::destroy() const
{
	delete _ref;
	delete[] _data;
}


template <class T>
inline uint FixedImage<T>::get_index(uint x, uint y, uint channel) const
{
	//return _size_x * (_size_y * channel + y) + x;	// NOTE: channel by channel
	return _number_of_channels * (_size_x * y + x) + channel;
}


/* ==================== Image ==================== */


template <class T>
Image<T>::Image()
 : FixedImage<T>()
{
}


template <class T>
Image<T>::Image(uint size_x, uint size_y)
 : FixedImage<T>(size_x, size_y)
{
}


template <class T>
Image<T>::Image(uint size_x, uint size_y, uint number_of_channels)
 : FixedImage<T>(size_x, size_y, number_of_channels)
{
}


template <class T>
Image<T>::Image(uint size_x, uint size_y, T default_value)
 : FixedImage<T>(size_x, size_y, default_value)
{
}


template <class T>
Image<T>::Image(uint size_x, uint size_y, uint number_of_channels, T default_value)
 : FixedImage<T>(size_x, size_y, number_of_channels, default_value)
{
}


template <class T>
Image<T>::Image(Shape size)
 : FixedImage<T>(size)
{
}


template <class T>
Image<T>::Image(Shape size, uint number_of_channels)
 : FixedImage<T>(size, number_of_channels)
{
}


template <class T>
Image<T>::Image(Shape size, T default_value)
 : FixedImage<T>(size, default_value)
{
}


template <class T>
Image<T>::Image(Shape size, uint number_of_channels, T default_value)
 : FixedImage<T>(size, number_of_channels, default_value)
{
}


template <class T>
Image<T>::Image(const Image<T> &source)
 : FixedImage<T>(source)	// NOTE: no data copying, ref++
{
}


template <class T>
Image<T>::Image(const FixedImage<T> &source)
{
	if (source._ref) {
		this->_size_x = source._size_x;
		this->_size_y = source._size_y;
		this->_number_of_channels = source._number_of_channels;
		FixedImage<T>::init(source._size_x, source._size_y, source._number_of_channels);
		memcpy(this->_data, source._data,  source._number_of_channels * source._size_y * source._size_x * sizeof(T));
	} else {
		this->_size_x = 0;
		this->_size_y = 0;
		this->_number_of_channels = 0;
		this->_ref = 0;
		this->_data = 0;
	}
}


template <class T>
Image<T>::~Image()
{
}


template <class T>
Image<T>& Image<T>::operator= (const Image<T> &other)
{
	// check for self-assignment
	if(this == &other) {
		return *this;
	}

	// finish all deals with the previous data
	Image<T>::release();

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

	return *this;
}


template <class T>
Image<T>& Image<T>::operator= (const FixedImage<T> &other)
{
	// finish all deals with the previous data
	Image<T>::release();

	if (other._ref) {
		this->_size_x = other._size_x;
		this->_size_y = other._size_y;
		this->_number_of_channels = other._number_of_channels;
		Image<T>::init(other._size_x, other._size_y, other._number_of_channels);
		memcpy(this->_data, other._data,  other._number_of_channels * other._size_y * other._size_x * sizeof(T));
	} else {
		this->_size_x = 0;
		this->_size_y = 0;
		this->_number_of_channels = 0;
		this->_ref = 0;
		this->_data = 0;
	}

	return *this;
}


template <class T>
bool Image<T>::operator== (const FixedImage<T> &other) const
{
	return this->_ref == other._ref;
}


template <class T>
bool Image<T>::operator== (const Image<T> &other) const
{
	return this->_ref == other._ref;
}


/**
 * Returns a reference to the element without range checking.
 */
template <class T>
T& Image<T>::operator() (uint x, uint y)
{
	return this->_data[Image<T>::get_index(x, y, 0)];
}


/**
 * Returns a reference to the element without range checking.
 */
template <class T>
T& Image<T>::operator() (uint x, uint y, uint channel)
{
	return this->_data[Image<T>::get_index(x, y, channel)];
}


/**
 * Returns a reference to the element without range checking.
 */
template <class T>
T& Image<T>::operator() (Point p)
{
	return this->_data[Image<T>::get_index(p.x, p.y, 0)];
}


/**
 * Returns a reference to the element without range checking.
 */
template <class T>
T& Image<T>::operator() (Point p, uint channel)
{
	return this->_data[Image<T>::get_index(p.x, p.y, channel)];
}


/**
 * Returns a reference to the element with range checking.
 * @note Throws std::out_of_range exception, if out of range.
 */
template <class T>
T& Image<T>::at(uint x, uint y)
{
	if (x >= this->_size_x || y >= this->_size_y) {
		throw std::out_of_range("x or y coordinate is out of range");
	}

	return this->_data[Image<T>::get_index(x, y, 0)];
}


/**
 * Returns a reference to the element with range checking.
 * @note Throws std::out_of_range exception, if out of range.
 */
template <class T>
T& Image<T>::at(uint x, uint y, uint channel)
{
	if (x >= this->_size_x || y >= this->_size_y || channel >= this->_number_of_channels) {
		throw std::out_of_range("channel, x or y coordinate is out of range");
	}

	return this->_data[Image<T>::get_index(x, y, channel)];
}


/**
 * Returns a reference to the element with range checking.
 * @note Throws std::out_of_range exception, if out of range.
 */
template <class T>
T& Image<T>::at(Point p)
{
	if (p.x < 0 || p.x >= this->_size_x || p.y < 0 || p.y >= this->_size_y) {
		throw std::out_of_range("x or y coordinate is out of range");
	}

	return this->_data[Image<T>::get_index(p.x, p.y, 0)];
}


/**
 * Returns a reference to the element with range checking.
 * @note Throws std::out_of_range exception, if out of range.
 */
template <class T>
T& Image<T>::at(Point p, uint channel)
{
	if (p.x < 0 || p.x >= this->_size_x || p.y < 0 || p.y >= this->_size_y || channel >= this->_number_of_channels) {
		throw std::out_of_range("channel, x or y coordinate is out of range");
	}

	return this->_data[Image<T>::get_index(p.x, p.y, channel)];
}


/**
 * Assigns a given value to all elements.
 */
template <class T>
void Image<T>::fill(const T &value)
{
	Image<T>::fill_internal(value);
}


/**
 * Returns pointer to internal data.
 */
template <class T>
T* Image<T>::raw()
{
	return this->_data;
}


/**
 * Invokes deep copy.
 */
template <class T>
Image<T> Image<T>::clone() const
{
	Image<T> clone;
	clone._size_x = this->_size_x;
	clone._size_y = this->_size_y;
	clone._number_of_channels = this->_number_of_channels;

	if (this->_ref) {
		clone.init(this->_size_x, this->_size_y, this->_number_of_channels);
		memcpy(clone._data, this->_data,  this->_number_of_channels * this->_size_y * this->_size_x * sizeof(T));
	}

	return clone;
}

/* Protected */

template <class T>
void Image<T>::destroy() const
{
	// no additional data
	FixedImage<T>::destroy();
}
