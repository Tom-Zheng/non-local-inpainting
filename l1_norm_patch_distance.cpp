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

#include "l1_norm_patch_distance.h"

L1NormPatchDistance::L1NormPatchDistance()
	: APatchDistance() { }

L1NormPatchDistance::L1NormPatchDistance(Shape &patch_size, float gaussian_sigma)
	: APatchDistance(patch_size, gaussian_sigma) { }


float L1NormPatchDistance::calculate(const Point &source_point,
							   	     const Point &target_point)
{
	if (_patch_weighting.is_empty()) {
		_patch_weighting = GaussianWeights::calculate(_patch_size.size_x,
													  _patch_size.size_y,
													  _gaussian_sigma,
													  _gaussian_sigma);
	}

	// calculate distance between source and target patches
	int radius_x = _patch_size.size_x / 2;
	int radius_y = _patch_size.size_y / 2;

	// NOTE: direct access - we sacrifice readability in favor of performance
	const float *source_points = _source.raw();
	const float *target_points = _target.raw();
	const float *p_weight = _patch_weighting.raw();
	int number_of_channels = _source.get_number_of_channels();
	int source_stride = _source.get_size_x();
	int target_stride = _target.get_size_x();


	float distance = 0.0;
	if (number_of_channels == 3) {
		// NOTE: unrolled channel loop for 3 channels
		for (int dy = -radius_y; dy < ((int)_patch_size.size_y - radius_y); dy++) {
			for (int dx = -radius_x; dx < ((int)_patch_size.size_x - radius_x); dx++, ++p_weight) {
				// NOTE: this is the same as: _patch_weighting(dx + radius_x, dy + radius_y);
				float weight = *p_weight;

				// calculate partial result as a weighted norm
				int source_id = 3 * (source_stride * (source_point.y + dy) + (source_point.x + dx));
				int target_id = 3 * (target_stride * (target_point.y + dy) + (target_point.x + dx));

				float norm = abs(source_points[source_id] - target_points[target_id]) +
							 abs(source_points[source_id + 1] - target_points[target_id + 1]) +
							 abs(source_points[source_id + 2] - target_points[target_id + 2]);

				distance += weight * norm;
			}
		}
	} else if (number_of_channels == 1) {
		// NOTE: no channel loop for 1 channel
		for (int dy = -radius_y; dy < ((int)_patch_size.size_y - radius_y); dy++) {
			for (int dx = -radius_x; dx < ((int)_patch_size.size_x - radius_x); dx++, ++p_weight) {
				// NOTE: this is the same as: _patch_weighting(dx + radius_x, dy + radius_y);
				float weight = *p_weight;

				// calculate partial result as a weighted norm
				int source_id = source_stride * (source_point.y + dy) + (source_point.x + dx);
				int target_id = target_stride * (target_point.y + dy) + (target_point.x + dx);

				float norm = abs(source_points[source_id] - target_points[target_id]);

				distance += weight * norm;
			}
		}
	} else {
		for (int dy = -radius_y; dy < ((int)_patch_size.size_y - radius_y); dy++) {
			for (int dx = -radius_x; dx < ((int)_patch_size.size_x - radius_x); dx++, ++p_weight) {
				// NOTE: this is the same as: _patch_weighting(dx + radius_x, dy + radius_y);
				float weight = *p_weight;

				// calculate partial result as a weighted norm
				float norm = 0.0;
				int source_id = number_of_channels * (source_stride * (source_point.y + dy) + (source_point.x + dx));
				int target_id = number_of_channels * (target_stride * (target_point.y + dy) + (target_point.x + dx));
				for (int ch = 0; ch < number_of_channels; ch++) {
					// NOTE: this is the same as: source(source_point + dp);
					float color_a = source_points[source_id + ch];
					// NOTE: this is the same as: target(target_point + dp);
					float color_b = target_points[target_id + ch];

					norm += abs(color_a - color_b);
				}

				distance += weight * norm;
			}
		}
	}

	return distance / (_patch_size.size_x * _patch_size.size_y);
}
