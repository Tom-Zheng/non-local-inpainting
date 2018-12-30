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

#include "distance_transform.h"

namespace DistanceTransform {

/**
 * For each point assigns the distance to the nearest point NOT belonging to the mask
 * (masked points receive non-zero distances).
 */
Image<float> calculate(FixedMask mask)
{
	// NOTE: it's hard to give meaningful names to the variables,
	// therefore, most of them are named like in the paper.

	// allocate buffer to store partial results ('G' function's values in the paper)
	Image<float> g(mask.get_size_x(), mask.get_size_y());

	// allocate memory for result
	Image<float> distances(mask.get_size_x(), mask.get_size_y());

	// NOTE: infinity is replaced by the maximum possible distance
	const float INF = mask.get_size_x() + mask.get_size_y();

	// first phase
	for (uint x = 0; x < mask.get_size_x(); x++) {
		// scan 1 (top to bottom)
		if (mask.get(x, 0)) {
			g(x, 0) = INF;
		} else {
			g(x, 0) = 0.0;
		}

		for (uint y = 1; y < mask.get_size_y(); y++) {
			if (mask.get(x, y)) {
				g(x, y) = g(x, y - 1) + 1.0;
			} else {
				g(x, y) = 0.0;
			}
		}

		// scan 2 (bottom to top)
		for (int y = mask.get_size_y() - 2; y >= 0; y--) {
			if (g(x, y + 1) < g(x, y)) {
				g(x, y) = g(x, y + 1) + 1.0;
			}
		}
	}

	// second phase
	vector<int> s(mask.get_size_x());
	vector<int> t(mask.get_size_x());
	for (uint y = 0; y < mask.get_size_y(); y++) {
		int q = 0;
		s[0] = 0;
		t[0] = 0;

		// scan 3 (left to right)
		for (uint u = 1; u < mask.get_size_x(); u++) {
			while (q >= 0 && _Details::f(t[q], s[q], g(s[q], y)) > _Details::f(t[q], u, g(u, y)) ) {
				q--;
			}

			if (q < 0) {
				q = 0;
				s[0] = u;
			} else {
				uint w = 1 + _Details::sep(s[q], u, g(s[q], y), g(u, y), INF);
				if (w < mask.get_size_x()) {
					q++;
					s[q] = u;
					t[q] = w;
				}
			}
		}

		// scan 4 (right to left)
		for (int u = mask.get_size_x() - 1; u >= 0; u--) {
			distances(u, y) = _Details::f(u, s[q], g(s[q], y));
			if (u == t[q]) {
				q--;
			}
		}

		s.clear();
		t.clear();
	}

	return distances;
}



namespace _Details {

inline float f(int x, int x_i, float g_i)
{
	return sqrt(pow((float)(x - x_i), 2) + pow(g_i, 2));
}


inline int sep(int i, int u, float g_i, float g_u, int inf)
{
	return (pow((float)u, 2) - pow((float)i, 2) + pow(g_u, 2) - pow(g_i, 2)) / (2 * (u - i));
}

} // namespace _Details

} // namespace DistanceTransform


