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

#ifndef L2_RGBD_PATCH_DISTANCE_H_
#define L2_RGBD_PATCH_DISTANCE_H_

#include "a_patch_distance.h"
#include "gradient.h"

/**
 * Implements L2-Norm patch distance calculation method
 * which takes into account not only intensity differences,
 * but also gradient differences. The impact of the intensity
 * term over the gradient term is controlled by the 'lambda'
 * parameter (lambda equal to 1.0 implies simple L2-Norm
 * with no gradient impact).
 */
class L2RGBDPatchDistance : public APatchDistance
{
public:
	L2RGBDPatchDistance();
	L2RGBDPatchDistance(float lambda_rgb, float lambda_d, Shape &patch_size, float gaussian_sigma);

	virtual ~L2RGBDPatchDistance() {};

	virtual void initialize(FixedImage<float> source, FixedImage<float> target);

	virtual float calculate(const Point &source_point,
							const Point &target_point);

private:
	FixedImage<float> _source_gradient;
	FixedImage<float> _target_gradient;
	float _lambda_rgb;
	float _lambda_d;
};




#endif /* L2_RGBD_PATCH_DISTANCE_H_ */
