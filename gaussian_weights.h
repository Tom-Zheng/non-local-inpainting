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

#ifndef GAUSSIAN_WEIGHTS_H_
#define GAUSSIAN_WEIGHTS_H_

#include "image.h"

/**
 * Contains methods for calculating discrete Gaussian kernel.
 */
namespace GaussianWeights {

static const unsigned int SEGMENTS = 1;	// NOTE: 1 means no integration

Image<float> calculate(unsigned int size_x, unsigned int size_y, float sigma_x, float sigma_y);
Image<float> calculate_2d(unsigned int size, float sigma);
Image<float> calculate_1d(unsigned int size, float sigma);
}	// namespace GaussianWeights



#endif /* GAUSSIAN_WEIGHTS_H_ */
