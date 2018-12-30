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

#include "gradient.h"

namespace Gradient {

/**
 * Calculates forward or backward gradient for a given image.
 */
Image<float> calculate(FixedImage<float> image, GradientType x_type, GradientType y_type)
{
	int size_x = image.get_size_x();
	int size_y = image.get_size_y();
	int number_of_channels = image.get_number_of_channels();

	// allocate memory
	Image<float> gradient(size_x, size_y, (uint)(number_of_channels * 2));
	float *gradient_data = gradient.raw();
	const float *image_data = image.raw();

	// calculate gradient
	for (int x = 0; x < size_x; x++) {
		for (int y = 0; y < size_y; y++) {
			int index = number_of_channels * (size_x * y + x);
			for (int ch = 0; ch < number_of_channels; ch++) {
				if (x_type == Forward) {
					if (x == size_x - 1) {
						gradient_data[(index + ch) * 2] = 0.0;
					} else {
						float a = image_data[index + ch];
						float b = image_data[index + number_of_channels + ch];
						gradient_data[(index + ch) * 2] = b - a;
					}
				} else {
					if (x == 0) {
						gradient_data[(index + ch) * 2] = 0.0;
					} else {
						float a = image_data[index - number_of_channels + ch];
						float b = image_data[index + ch];
						gradient_data[(index + ch) * 2] = b - a;
					}
				}

				if (y_type == Forward) {
					if (y == size_y - 1) {
						gradient_data[(index + ch) * 2 + 1] = 0.0;
					} else {
						float a = image_data[index + ch];
						float b = image_data[index + number_of_channels * size_x + ch];
						gradient_data[(index + ch) * 2 + 1] = b - a;
					}
				} else {
					if (y == 0) {
						gradient_data[(index + ch) * 2 + 1] = 0.0;
					} else {
						float a = image_data[index - number_of_channels * size_x + ch];
						float b = image_data[index + ch];
						gradient_data[(index + ch) * 2 + 1] = b - a;
					}
				}
			}	// for (int ch = 0; ...)
		}	// for (int y = 0; ...)
	}	// for (int x = 0; ...)

	return gradient;
}

}
