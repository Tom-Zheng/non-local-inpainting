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

#ifndef IMAGE_H_
#define IMAGE_H_

#include <cstring>
#include <algorithm>
#include <stdexcept>
#include "point.h"
#include "shape.h"

typedef unsigned int uint;

using namespace std;

template <class T>
class Image;	// forward declaration

/**
 * Container for a 2d image with point type T. Manages memory internally by
 * references counting. Constructing from the same type and assignment of
 * a value of the same type lead to a data sharing. Method clone() should be
 * used for explicit deep copy invocation.
 *
 * @note By design class provides no capabilities for changing its data,
 * therefore, 'Fixed' here should be considered as 'Immutable'.
 */
template <class T = float>
class FixedImage
{
friend class Image<T>;
public:
	FixedImage();
	FixedImage(uint size_x, uint size_y);
	FixedImage(uint size_x, uint size_y, uint number_of_channels);
	FixedImage(uint size_x, uint size_y, T default_value);
	FixedImage(uint size_x, uint size_y, uint number_of_channels, T default_value);
	FixedImage(Shape size);
	FixedImage(Shape size, uint number_of_channels);
	FixedImage(Shape size, T default_value);
	FixedImage(Shape size, uint number_of_channels, T default_value);
	FixedImage(const FixedImage<T> &source);		// without data copying, ref++
	FixedImage(const Image<T> &source);				// without data copying, ref++
	virtual ~FixedImage();

	FixedImage<T>& operator= (const FixedImage<T> &other);		// without data copying, ref++
	FixedImage<T>& operator= (const Image<T> &other);			// without data copying, ref++

	bool operator== (const FixedImage<T> &other) const;
	bool operator== (const Image<T> &other) const;

	/// Is current image [not] empty.
	bool is_empty() const;
	bool is_not_empty() const;

	uint get_size_x() const;
	uint get_size_y() const;
	Shape get_size() const;
	uint get_number_of_channels() const;

	/// Returns read-only value without range checking.
	const T& operator() (uint x, uint y) const;
	const T& operator() (uint x, uint y, uint channel) const;
	const T& operator() (Point p) const;
	const T& operator() (Point p, uint channel) const;

	/// Returns read-only value with range checking.
	/// Throws std::out_of_range exception, if out of range.
	const T& at(uint x, uint y) const;
	const T& at(uint x, uint y, uint channel) const;
	const T& at(Point p) const;
	const T& at(Point p, uint channel) const;

	/// Checks if all indexes are in range and modifies the value parameter.
	bool try_get_value(uint x, uint y, T& value) const;
	bool try_get_value(uint x, uint y, uint channel, T& value) const;
	bool try_get_value(Point p, T& value) const;
	bool try_get_value(Point p, uint channel, T& value) const;

	/// Returns pointer to internal data.
	const T* raw() const;

	/// Invokes deep copy.
	FixedImage<T> clone() const;

protected:
	struct __Ref {
		int counter;
	};

	uint _size_x, _size_y;
	uint _number_of_channels;
	T* _data;
	__Ref *_ref;

	inline void init(uint size_x, uint size_y, uint number_of_channels);
	void fill_internal(const T &value);

	void release() const;
	virtual void destroy() const;

	inline uint get_index(uint x, uint y, uint channel) const;
};


/**
 * Container for a 2d image with point type T. Manages memory internally by
 * references counting. Constructing from the same type and assignment of
 * a value of the same type lead to a data sharing. Method clone() should be
 * used for explicit deep copy invocation.
 *
 * @note Extends the FixedImage<T> class with the data modification capabilities.
 */
template <class T>
class Image : public FixedImage<T>
{
friend class FixedImage<T>;
public:
	Image();
	Image(uint size_x, uint size_y);
	Image(uint size_x, uint size_y, uint number_of_channels);
	Image(uint size_x, uint size_y, T default_value);
	Image(uint size_x, uint size_y, uint number_of_channels, T default_value);
	Image(Shape size);
	Image(Shape size, uint number_of_channels);
	Image(Shape size, T default_value);
	Image(Shape size, uint number_of_channels, T default_value);
	Image(const Image<T> &source);			// without data copying, ref++
	Image(const FixedImage<T> &source);		// deep copy
	virtual ~Image();

	Image<T>& operator= (const Image<T> &other);			// without data copying, ref++
	Image<T>& operator= (const FixedImage<T> &other);		// deep copy

	bool operator== (const FixedImage<T> &other) const;
	bool operator== (const Image<T> &other) const;

	/// Returns a reference to the element without range checking.
	using FixedImage<T>::operator();
	T& operator() (uint x, uint y);
	T& operator() (uint x, uint y, uint channel);
	T& operator() (Point p);
	T& operator() (Point p, uint channel);

	/// Returns a reference to the element with range checking.
	/// Throws std::out_of_range exception, if out of range.
	using FixedImage<T>::at;
	T& at(uint x, uint y);
	T& at(uint x, uint y, uint channel);
	T& at(Point p);
	T& at(Point p, uint channel);

	/// Assigns a given value to all elements.
	void fill(const T &value);

	/// Returns pointer to internal data.
	T* raw();

	/// Invokes deep copy.
	Image<T> clone() const;

protected:
	virtual void destroy() const;
};

// NOTE: include implementation, because Image is a template
#include "image.hpp"

#endif /* IMAGE_H_ */
