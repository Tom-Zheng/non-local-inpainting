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

#include "patch_non_local_medians.h"

PatchNonLocalMedians::PatchNonLocalMedians()
: AImageUpdating() { }


PatchNonLocalMedians::PatchNonLocalMedians(Shape patch_size, float gaussian_sigma)
: AImageUpdating(patch_size, gaussian_sigma) { }


double PatchNonLocalMedians::update(Image<float> image,
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

		// define borders of the patch
		int x_a = max(x - half_patch_size_x, 0);
		int x_b = min(x + half_patch_size_x, (int)image_size.size_x - 1);
		int y_a = max(y - half_patch_size_y, 0);
		int y_b = min(y + half_patch_size_y, (int)image_size.size_y - 1);

		float total_weight = 0.0;
				
		// collect all candidate colors with their weights
		vector< vector< pair< float, float> > > color_weights_map(number_of_channels);
		for (int i=0; i<number_of_channels; i++) 
			color_weights_map[i].reserve(_patch_size.size_x * _patch_size.size_y);

		for (int i = x_a; i <= x_b; i++) {
			for (int j = y_a; j <= y_b; j++) {
				// get center of the corresponding patch
				Point contributor = nnf(i, j);

				if (contributor.x == -1) {
					continue;
				}

				// apply offset w.r.t the patch center to obtain contributing point
				contributor.x += x - i;
				contributor.y += y - j;

				if (image_size.contains(contributor)) {
					// get weight for the contributor
					float weight = _patch_weighting(i - x_a, j - y_a);
					if (confidence_mask.is_not_empty() && inpainting_domain.test(i, j)) {	// NOTE: outside the inpainting domain Confidence is 1.0, therefore multiplication might be skipped
						weight *= confidence_mask(i, j);
					}
					total_weight += weight;

					for (int ch = 0; ch < number_of_channels; ch++) {
						color_weights_map[ch].push_back(  pair< float,float> (image(contributor, ch), weight) ); 
					}
				}
			}
		}

		if (color_weights_map[0].size() == 0)
			continue;

		// NOTE: three channels are processed separately
		for (int ch = 0; ch < number_of_channels; ch++) {
			float color_value =	compute_weighted_median(color_weights_map[ch], total_weight);

			// add to total difference
			float prev_value = image(x, y, ch);
			total_difference += (prev_value - color_value) * (prev_value - color_value);

			image(x, y, ch) = color_value;
		}

	}

	return total_difference;
}

/* Private */

float PatchNonLocalMedians::compute_weighted_median(vector< pair< float, float> > &weight_map, float total_weight)
{
	std::sort(weight_map.begin(), weight_map.end(), CompareByColor());
	
	float median_color = 0.0;
	float half_total_weight = total_weight / 2;
	vector<pair<float, float> >::iterator it;
	for (it = weight_map.begin(); it != weight_map.end(); ++it) {
		total_weight -= it->second;
		if (total_weight < half_total_weight) {
			median_color = it->first;
			break;
		}
	}

	return median_color;
}


