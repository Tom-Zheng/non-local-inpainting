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

#ifndef L2_NORM_PATCH_DISTANCE_H_
#define L2_NORM_PATCH_DISTANCE_H_

#include <map>
#include "a_patch_distance.h"

using namespace std;

/**
 * Implements L2-Norm patch distance calculation method.
 */
class L2NormPatchDistance : public APatchDistance
{
public:
	L2NormPatchDistance();
	L2NormPatchDistance(Shape &patch_size, float gaussian_sigma);

	virtual ~L2NormPatchDistance() {};

	virtual float calculate(const Point &source_point,
							const Point &target_point);
};


#endif /* L2_NORM_PATCH_DISTANCE_H_ */
