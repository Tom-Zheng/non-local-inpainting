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

#ifndef PATCH_NON_LOCAL_MEANS_H_
#define PATCH_NON_LOCAL_MEANS_H_

#include "a_image_updating.h"
#include "image.h"

using namespace std;

/**
 * Implements Non-Local Means image updating scheme.
 */
class PatchNonLocalMeans : public AImageUpdating
{
public:
	PatchNonLocalMeans();
	PatchNonLocalMeans(Shape patch_size, float gaussian_sigma);
	virtual ~PatchNonLocalMeans() {}

	virtual double update(Image<float> image,
						  Image<float> original_image,
						  FixedMask inpainting_domain,
						  FixedMask extended_inpainting_domain,
						  FixedImage<Point> nnf,
						  FixedImage<float> confidence_mask);

};


#endif /* PATCH_NON_LOCAL_MEANS_H_ */
