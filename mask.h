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

#ifndef MASK_H_
#define MASK_H_

#include <vector>
#include <stdexcept>

#include "image.h"
#include "mask_iterator.h"

using namespace std;

/// forward declaration
class Mask;
class MaskIterator;

/**
 * Container for a 2d binary mask based on the FixedImage<bool> class.
 * Provides capabilities for iterating through the masked points and
 * retrieving masked points as a vector.
 *
 * @note By design class provides no capabilities for changing its data,
 * therefore, 'Fixed' here should be considered as 'Immutable'.
 */
class FixedMask : public FixedImage<bool>
{
friend class Mask;
public:
	typedef MaskIterator iterator;

	FixedMask();
	FixedMask(uint size_x, uint size_y);
	FixedMask(uint size_x, uint size_y, bool default_value);
	FixedMask(Shape size);
	FixedMask(Shape size, bool default_value);
	FixedMask(const Image<bool> &source);			// deep copy
	FixedMask(const FixedImage<bool> &source);		// deep copy
	FixedMask(const Mask &source);					// without data copying, ref++
	FixedMask(const FixedMask &source);				// without data copying, ref++
	virtual ~FixedMask();

	FixedMask& operator= (const Mask &other);			// without data copying, ref++
	FixedMask& operator= (const FixedMask &other);		// without data copying, ref++

	iterator begin() const;
	iterator end() const;
	iterator rbegin() const;
	iterator rend() const;

	/// Returns the element without range checking. Does not affect internal cache.
	bool get(uint x, uint y) const;
	bool get(uint x, uint y, uint channel) const;
	bool get(Point p)  const;
	bool get(Point p, uint channel)  const;

	/// Returns the element with range checking (returns 'false', if out of range). Does not affect internal cache.
	bool test(uint x, uint y)  const;
	bool test(uint x, uint y, uint channel)  const;
	bool test(Point p)  const;
	bool test(Point p, uint channel)  const;

	/// Invokes deep copy.
	FixedMask clone() const;

	/// Invokes deep copy and then invert the cloned mask.
	FixedMask clone_invert() const;

	/// Returns masked points as a vector.
    vector<Point> get_masked_points() const;

	/// Methods used by MaskIterator
	Point first() const;
	Point last() const;
	Point next(const Point &current) const;
	Point prev(const Point &current) const;

	/// Bounding box
	Point bounding_box_top_left() const;
	Point bounding_box_bottom_right() const;

protected:
	struct __Internal
	{
		Point first;
		Point last;
		Point top_left;			// bounding box
		Point bottom_right;		// bounding box
		bool is_first_last_valid;
		bool is_bounding_box_valid;
		bool is_points_cache_valid;
		vector<Point> points_cache;
	};

	mutable __Internal *_internal;

	virtual void destroy() const;

	inline void init_internal(Point first, Point last, bool is_first_last_valid, Point top_left, Point bottom_rigth, bool is_bounding_box_valid);
	inline void actualize_first_last() const;
	inline void actualize_points_cache() const;
};


/**
 * Container for a 2d binary mask based on the FixedImage<bool> class.
 * Provides capabilities for iterating through the masked points and
 * retrieving masked points as a vector.
 *
 * @note Extends the FixedMask class with the data modification capabilities.
 */
class Mask : public FixedMask
{
friend class FixedMask;
public:
	typedef MaskIterator iterator;

	Mask();
	Mask(uint size_x, uint size_y);
	Mask(uint size_x, uint size_y, bool default_value);
	Mask(Shape size);
	Mask(Shape size, bool default_value);
	Mask(const Image<bool> &source);			// deep copy
	Mask(const FixedImage<bool> &source);		// deep copy
	Mask(const Mask &source);					// without data copying, ref++
	Mask(const FixedMask &source);				// deep copy
	virtual ~Mask();

	Mask& operator= (const Mask &other);		// without data copying, ref++
	Mask& operator= (const FixedMask &other);	// deep copy

	/// Returns a reference to the element without range checking. Invalidates internal cache.
	using FixedImage<bool>::operator();
	bool& operator() (uint x, uint y);
	bool& operator() (uint x, uint y, uint channel);
	bool& operator() (Point p);
	bool& operator() (Point p, uint channel);

	/// Returns a reference to the element with range checking.
	/// Throws std::out_of_range exception, if out of range. Invalidates internal cache.
	using FixedImage<bool>::at;
	bool& at(uint x, uint y);
	bool& at(uint x, uint y, uint channel);
	bool& at(Point p);
	bool& at(Point p, uint channel);

	/// Sets element to 'true' with range checking (does nothing, if out of range).
	/// Invalidates internal cache.
	void mask(uint x, uint y);
	void mask(uint x, uint y, uint channel);
	void mask(Point p);
	void mask(Point p, uint channel);

	/// Sets element to 'false' with range checking (does nothing, if out of range).
	/// Invalidates internal cache.
	void unmask(uint x, uint y);
	void unmask(uint x, uint y, uint channel);
	void unmask(Point p);
	void unmask(Point p, uint channel);

	/// Invokes deep copy.
	Mask clone() const;

	/// Invokes deep copy and then invert the cloned mask.
	Mask clone_invert() const;

	/// Invert current mask.
    void invert();

protected:

    virtual void destroy() const;
};

#endif /* MASK_H_ */
