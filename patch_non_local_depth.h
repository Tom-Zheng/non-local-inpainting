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

#ifndef PATCH_NON_LOCAL_DEPTH_H_
#define PATCH_NON_LOCAL_DEPTH_H_

#include <cstring>
#include "a_image_updating.h"
#include "patch_non_local_poisson.h"
#include "patch_non_local_means.h"
#include "io_utility.h"

using namespace std;

/**
 * Implements image updating scheme for RGBD inpainting
 */
class PatchNonLocalDepth : public AImageUpdating
{
public:
	PatchNonLocalDepth();

	PatchNonLocalDepth(Shape patch_size,
						 float gaussian_sigma,
						 float lambda_rgb,
						 float lambda_d,
						 float conjugate_gradient_tolerance = 1e-10,
						 int conjugate_gradient_iterations_limit = 1000);

	virtual ~PatchNonLocalDepth() {delete rgb_updating;    delete depth_updating;}

	virtual double update(Image<float> image,
						  Image<float> original_image,
						  FixedMask inpainting_domain,
						  FixedMask extended_inpainting_domain,
						  FixedImage<Point> nnf,
						   FixedImage<float> confidence_mask);

private:
    AImageUpdating *rgb_updating;
	AImageUpdating *depth_updating;
};




#endif /* PATCH_NON_LOCAL_DEPTH_H_ */
