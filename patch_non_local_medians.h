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

#ifndef PATCH_NON_LOCAL_MEDIANS_H_
#define PATCH_NON_LOCAL_MEDIANS_H_

#include <map>
#include <list>
#include "a_image_updating.h"

using namespace std;

/**
 * Implements Non-Local Medians image updating scheme.
 */
class PatchNonLocalMedians : public AImageUpdating
{
public:
	PatchNonLocalMedians();
	PatchNonLocalMedians(Shape patch_size, float gaussian_sigma);
	virtual ~PatchNonLocalMedians() {}

	virtual double update(Image<float> image,
						  Image<float> original_image,
						  FixedMask inpainting_domain,
						  FixedMask extended_inpainting_domain,
						  FixedImage<Point> nnf,
						  FixedImage<float> confidence_mask);

private:
	// function object for comparing two pairs by the first elements (which is color in our case)
	struct CompareByColor : public binary_function<pair<float, float>, pair<float, float>, bool>
	{
		inline bool operator()(const pair<float, float>& first, const pair<float, float>& second)
		{
			return first.first < second.first;
		}
	};

	float compute_weighted_median(vector< pair< float, float> > &weight_map, float total_weight);

};

#endif /* PATCH_NON_LOCAL_MEDIANS_H_ */
