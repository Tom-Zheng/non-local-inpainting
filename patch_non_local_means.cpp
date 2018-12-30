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

#include "patch_non_local_means.h"

PatchNonLocalMeans::PatchNonLocalMeans()
: AImageUpdating() { }


PatchNonLocalMeans::PatchNonLocalMeans(Shape patch_size, float gaussian_sigma)
: AImageUpdating(patch_size, gaussian_sigma) { }


double PatchNonLocalMeans::update(Image<float> image,
								  Image<float> original_image,
								  FixedMask inpainting_domain,
								  FixedMask extended_inpainting_domain,
								  FixedImage<Point> nnf,
								  FixedImage<float> confidence_mask)
{
	int half_patch_size_x = _patch_size.size_x / 2;
	int half_patch_size_y = _patch_size.size_y / 2;
	Shape image_size = image.get_size();
	int number_of_channels = image.get_number_of_channels();

	if (_patch_weighting.is_empty()) {
		_patch_weighting = GaussianWeights::calculate(_patch_size.size_x,
													  _patch_size.size_y,
													  _gaussian_sigma,
													  _gaussian_sigma);

	}

	double total_difference = 0.0;
	Mask::iterator it;
	for (it = inpainting_domain.begin(); it != inpainting_domain.end(); ++it) {
		int x = it->x;
		int y = it->y;

		int x_a = max(x - half_patch_size_x, 0);
		int x_b = min(x + half_patch_size_x, (int)image_size.size_x - 1);
		int y_a = max(y - half_patch_size_y, 0);
		int y_b = min(y + half_patch_size_y, (int)image_size.size_y - 1);

		float total_value[number_of_channels];
		for (int ch = 0; ch < number_of_channels; ch++) {
			total_value[ch] = 0.0;
		}

		float total_weight = 0.0;
		for (int i = x_a; i <= x_b; i++) {
			for (int j = y_a; j <= y_b; j++) {
				// get center of the corresponding patch
				Point contributor = nnf(i, j);

				if (contributor.x != -1) {
					// apply offset w.r.t the patch center to obtain contributing point
					contributor.x += x - i;
					contributor.y += y - j;

					if (image_size.contains(contributor)) {
						float weight = _patch_weighting(i - x_a, j - y_a);
						if (confidence_mask.is_not_empty() && inpainting_domain.test(i, j)) {	// NOTE: outside the inpainting domain Confidence is 1.0, therefore multiplication might be skipped
							weight *= confidence_mask(i, j);
						}
						total_weight += weight;

						for (int ch = 0; ch < number_of_channels; ch++) {
							total_value[ch] += image(contributor.x, contributor.y, ch) * weight;
						}
					}
				}
			}
		}

		if (total_weight > 0.0) {
			for (int ch = 0; ch < number_of_channels; ch++) {
				float color_value = total_value[ch] / total_weight;

				// add to the total difference
				float prev_value = image(x, y, ch);
				total_difference += (prev_value - color_value) * (prev_value - color_value);

				image(x, y, ch) = color_value;
			}
		}

	}

	return total_difference;
}
