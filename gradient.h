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

#ifndef GRADIENT_H_
#define GRADIENT_H_

#include "image.h"

namespace Gradient {

enum GradientType {
	Forward,
	Backward
};

Image<float> calculate(FixedImage<float> image, GradientType x_type = Forward, GradientType y_type = Forward);

}

#endif /* GRADIENT_H_ */
