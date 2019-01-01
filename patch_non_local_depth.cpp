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

#include "patch_non_local_depth.h"

PatchNonLocalDepth::PatchNonLocalDepth()
	: AImageUpdating()
{
	rgb_updating = new PatchNonLocalMeans();
	depth_updating = new PatchNonLocalPoisson();
}


PatchNonLocalDepth::PatchNonLocalDepth(Shape patch_size,
										   float gaussian_sigma,
										   float lambda,
										   float conjugate_gradient_tolerance,
										   int conjugate_gradient_iterations_limit)
	: AImageUpdating(patch_size, gaussian_sigma)
{
	// rgb_updating = new PatchNonLocalPoisson(patch_size, gaussian_sigma, lambda, 0.000001, 1000);
	rgb_updating = new PatchNonLocalMeans(patch_size, gaussian_sigma);
	depth_updating = new PatchNonLocalPoisson(patch_size, gaussian_sigma, lambda, 0.000001, 1000);
}

double PatchNonLocalDepth::update(Image<float> image,
									Image<float> original_image,
									FixedMask inpainting_domain,
									FixedMask extended_inpainting_domain,
									FixedImage<Point> nnf,
									FixedImage<float> confidence_mask)
{
	double total_difference = 0;
	Image<float> rgb(image.get_size(), (uint)3), rgb_orig(image.get_size(), (uint)3);
	Image<float> depth(image.get_size());
	Image<float> depth_orig(image.get_size());
	// Separate RGB and depth channels
	IOUtility::separate(image, rgb, depth);
	IOUtility::separate(image, rgb_orig, depth_orig);
	// Do updating accordingly
    total_difference += rgb_updating->update(image, original_image, inpainting_domain, extended_inpainting_domain, nnf, confidence_mask);
    total_difference += depth_updating->update(depth, depth_orig, inpainting_domain, extended_inpainting_domain, nnf, confidence_mask);
    // Merge channels
    image = IOUtility::cat(rgb, depth);
	return total_difference;
}
