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

#ifndef DISTANCE_TRANSFORM_H_
#define DISTANCE_TRANSFORM_H_

#include <vector>

#include "image.h"
#include "mask.h"

using namespace std;

/**
 * Contains implementation of "A General Algorithm for Computing
 * Distance Transforms in Linear Time" by Meijster et al.
 */
namespace DistanceTransform {

/// For each point assigns the distance to the nearest point NOT belonging to the mask
Image<float> calculate(FixedMask mask);


// NOTE: The following internal namespace contains implementation details. Do not call its members directly.
namespace _Details {

inline float f(int x, int x_i, float g_i);
inline int sep(int i, int u, float g_i, float g_u, int inf);

} // namespace _Details

} // namespace DistanceTransform




#endif /* DISTANCE_TRANSFORM_H_ */
