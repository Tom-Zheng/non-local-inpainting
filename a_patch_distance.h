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

#ifndef A_PATCH_DISTANCE_H_
#define A_PATCH_DISTANCE_H_

#include "gaussian_weights.h"
#include "image.h"

/**
 * Abstract base class for different patch distance calculation
 * methods (e.g. L2-Norm, L1-Norm, etc.). Contains functionality
 * related to the intra-patch weighting (using gaussian weights).
 */
class APatchDistance
{
public:
	APatchDistance();
	APatchDistance(Shape &patch_size, float gaussian_sigma);

	virtual ~APatchDistance();

	virtual void initialize(FixedImage<float> source, FixedImage<float> target);

	virtual float calculate(const Point &source_point,
							const Point &target_point) = 0;

	/// getters and setters for parameters
	float get_gaussian_sigma();
	void set_gaussian_sigma(float gaussian_sigma);
	Shape get_patch_size();
	void set_patch_size(Shape patch_size);

protected:
	FixedImage<float> _source;
	FixedImage<float> _target;
	Image<float> _patch_weighting;
	Shape _patch_size;
	float _gaussian_sigma;
};


#endif /* A_PATCH_DISTANCE_H_ */
