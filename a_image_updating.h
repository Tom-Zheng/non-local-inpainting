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

#ifndef A_IMAGE_UPDATING_H_
#define A_IMAGE_UPDATING_H_

#include "gaussian_weights.h"
#include "image.h"
#include "mask.h"
#include "shape.h"

/**
 * Abstract base class for different image updating methods
 * (e.g. Non-Local Means, Non-Local Poisson, etc.). Contains
 * functionality related to the intra-patch weighting (using
 * gaussian weights).
 */
class AImageUpdating
{
public:
	AImageUpdating();
	AImageUpdating(Shape patch_size, float gaussian_sigma);
	virtual ~AImageUpdating();

	virtual double update(Image<float> image,
						  Image<float> orig_image,
						  FixedMask inpainting_domain,
						  FixedMask extended_inpainting_domain,
						  FixedImage<Point> nnf,
						  FixedImage<float> confidence_mask) = 0;

	/// getters and setters for parameters
	float get_gaussian_sigma();
	void set_gaussian_sigma(float gaussian_sigma);
	Shape get_patch_size();
	void set_patch_size(Shape patch_size);

protected:
	float _gaussian_sigma;
	Shape _patch_size;
	Image<float> _patch_weighting;
};




#endif /* A_IMAGE_UPDATING_H_ */
