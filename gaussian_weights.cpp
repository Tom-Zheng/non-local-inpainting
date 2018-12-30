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

#include "gaussian_weights.h"

namespace GaussianWeights {

inline int get_index(int x, int y, unsigned int size_x);

Image<float> calculate(unsigned int size_x, unsigned int size_y, float sigma_x, float sigma_y)
{
	Image<float> weights(size_x, size_y, 0.0f);
	float *values = weights.raw();

	double x_a = - (float)size_x / 2;
	double y_a = - (float)size_y / 2;
	double step = 1.0 / SEGMENTS;
	double half_step = 0.5 / SEGMENTS;

	// 'Integrate' first half in 1D case, left-top corner in 2D case, etc. (including the central element)
	if (size_y == 1) {
		// 1D case
		for (unsigned int i = 0; i <= size_x / 2; i++) {
			double value = 0.0;

			for (unsigned int s = 0; s < SEGMENTS; s++) {
				double x = x_a + i + step * s + half_step;
				double partial = exp(-(pow(x, 2) / pow(sigma_x, 2)));
				value += partial;
			}
			values[get_index(i, 0, size_x)] = value;
		}
	} else {
		// 2D case
		for (unsigned int i = 0; i <= size_x / 2; i++) {
			for (unsigned int j = 0; j <= size_y / 2; j++) {
				double value = 0.0;
				for (unsigned int s = 0; s < SEGMENTS; s++) {
					double x = x_a + i + step * s + half_step;
					double y = y_a + j + step * s + half_step;
					double partial = exp(-(pow(x, 2) / (2 * pow(sigma_x, 2)) + pow(y, 2) / (2 * pow(sigma_y, 2))));
					value += partial;
				}
				values[get_index(i, j, size_x)] = value;
			}
		}
	}

	// Propagate values to the rest parts
	for (unsigned int i = 0; i <= size_x / 2; i++) {
		for (unsigned int j = 0; j <= size_y / 2; j++) {
			float value = values[get_index(i, j, size_x)];
			values[get_index(size_x - 1 - i, j, size_x)] = value;
			if (size_y > 1) {
				values[get_index(i, size_y - 1 - j, size_x)] = value;
				values[get_index(size_x - 1 - i, size_y - 1 - j, size_x)] = value;
			}
		}
	}

	// Normalize
	float sum = 0.0;
	for (unsigned int i = 0; i < size_x * size_y; i++) {
		sum += values[i];
	}
	for (unsigned int i = 0; i < size_x * size_y; i++) {
		values[i] /= sum;
	}

	return weights;
}


Image<float> calculate_2d(unsigned int size, float sigma)
{
	return calculate(size, size, sigma, sigma);
}


Image<float> calculate_1d(unsigned int size, float sigma)
{
	return calculate(size, 1, sigma, 0);
}


inline int get_index(int x, int y, unsigned int size_x)
{
	return size_x * y + x;
}

}	// namespace GaussianWeights




