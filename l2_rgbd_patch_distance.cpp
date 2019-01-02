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

#include "l2_rgbd_patch_distance.h"

L2RGBDPatchDistance::L2RGBDPatchDistance()
: APatchDistance()
{
	_lambda_rgb = 0.5;
	_lambda_d = 0.5;
}


L2RGBDPatchDistance::L2RGBDPatchDistance(float lambda_rgb, float lambda_d, Shape &patch_size, float gaussian_sigma)
: APatchDistance(patch_size, gaussian_sigma)
{
	_lambda_rgb = lambda_rgb;
	_lambda_d = lambda_d;
}


void L2RGBDPatchDistance::initialize(FixedImage<float> source, FixedImage<float> target)
{
	// calculate gradients
	_source_gradient = Gradient::calculate(source);
	_target_gradient = Gradient::calculate(target);

	APatchDistance::initialize(source, target);
}


float L2RGBDPatchDistance::calculate(const Point &source_point,
										 const Point &target_point)
{
	// get raw gradient values
	const float *source_gradient_values = _source_gradient.raw();
	const float *target_gradient_values = _target_gradient.raw();

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

	double distance = 0.0;
	double distance_depth = 0.0;
	double gradient_distance = 0.0;
	double gradient_distance_depth = 0.0;
	if (number_of_channels == 4)  {
		for (int dy = -radius_y; dy < ((int)_patch_size.size_y - radius_y); dy++) {
			for (int dx = -radius_x; dx < ((int)_patch_size.size_x - radius_x); dx++, ++p_weight) {
				// NOTE: this is the same as: _patch_weighting(dx + radius_x, dy + radius_y);
				float weight = *p_weight;

				// calculate partial result as a weighted norm
				double norm = 0.0;
				double norm_depth = 0.0;
				double gradient_norm = 0.0;
				double gradient_norm_depth = 0.0;
				int source_id = number_of_channels * (source_stride * (source_point.y + dy) + (source_point.x + dx));
				int target_id = number_of_channels * (target_stride * (target_point.y + dy) + (target_point.x + dx));
				for (int ch = 0; ch < number_of_channels; ch++) {
					// NOTE: this is the same as: source(source_point + dp);
					float color_a = source_points[source_id + ch];
					// NOTE: this is the same as: target(target_point + dp);
					float color_b = target_points[target_id + ch];

					float gradient_x_a = source_gradient_values[(source_id + ch) * 2];
					float gradient_x_b = target_gradient_values[(target_id + ch) * 2];
					float gradient_y_a = source_gradient_values[(source_id + ch) * 2 + 1];
					float gradient_y_b = target_gradient_values[(target_id + ch) * 2 + 1];
					
					if (ch == 3) {
						norm_depth += (color_a - color_b) * (color_a - color_b);
						gradient_norm_depth += (gradient_x_a - gradient_x_b) * (gradient_x_a - gradient_x_b) +
										 (gradient_y_a - gradient_y_b) * (gradient_y_a - gradient_y_b);
					} else {
						norm += (color_a - color_b) * (color_a - color_b);
						gradient_norm += (gradient_x_a - gradient_x_b) * (gradient_x_a - gradient_x_b) +
										 (gradient_y_a - gradient_y_b) * (gradient_y_a - gradient_y_b);
					}
				}
				distance += weight * norm;
				distance_depth += weight * norm_depth;
				gradient_distance += weight * gradient_norm;
				gradient_distance_depth += weight * gradient_norm_depth;
			}
		}
	} else {
		throw std::runtime_error("L2 RGBD Distance: Channels do not match.");
	}
	return (_lambda_rgb * distance + (1 - _lambda_rgb) * gradient_distance
	      + _lambda_d * distance_depth + (1 - _lambda_d) * gradient_distance_depth) 
		  / (_patch_size.size_x * _patch_size.size_y);
}
