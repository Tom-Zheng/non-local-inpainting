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

#ifndef SAMPLING_H_
#define SAMPLING_H_

#include <algorithm>
#include <math.h>
#include <utility> // for pair

#include "image.h"
#include "mask.h"
#include "gaussian_weights.h"
#include "shape.h"

using namespace std;

typedef unsigned int uint;

/**
 * Contains methods for downsampling both images and masks and
 * upsampling images.
 */
namespace Sampling {

// Downsample image by specified factor. Image is first filtered with
// appropriate Gaussian filtered. Bilinear interpolation is used to subsample
// filtered image.
Image<float> downsample(FixedImage<float> in, float factor);

// Downsample image and mask jointly by specified factor. The mask affects
// how the image is smoothed: the image is smoothed only where the mask is 1.
// Pixels out of the mask are treated as if they were out of the image. o
// specified sample_size. Image is multiplied by mask and filtered with
// appropriate Gaussian. Mask is also filtered. 
//
// The output image is defined by dividing the filtered masked image by the
// filtered mask (ie. the convolution is normalized by the number of pixels
// in the mask). 
//
// Bilinear interpolation is used to subsample filtered images.
std::pair<Image<float>, Mask> downsample_with_mask(FixedImage<float> image_in,
		Mask mask_in, float factor, float threshold = 0.5, bool mask_val = true);

// Downsample mask by specified factor. Mask is first filtered with
// appropriate Gaussian filtered. Bilinear interpolation is used to subsample
// filtered image. The resulting image is thresholded to obtain a 
// binary mask.
Mask downsample(FixedMask in, float factor, float threshold = 0.0);

Image<float> upsample(FixedImage<float> in, Shape size);


// NOTE: The following internal namespace contains implementation details. Do not call its members directly.
namespace _Details {

// Downsample image to specified sample_size. Image is first filtered with
// appropriate Gaussian filtered. Bilinear interpolation is used to subsample
// filtered image.
void downsample_internal(const float* in, float* out, uint size_x, uint size_y, uint sample_size_x, uint sample_size_y);

// Downsample image and mask jointly to specified sample_size. The mask affects
// how the image is smoothed: the image is smoothed only where the mask is 1.
// Pixels out of the mask are treated as if they were out of the image. o
// specified sample_size. Image is multiplied by mask and filtered with
// appropriate Gaussian. Mask is also filtered. 
//
// The output image is defined by dividing the filtered masked image by the
// filtered mask (ie. the convolution is normalized by the number of pixels
// in the mask). 
//
// Bilinear interpolation is used to subsample filtered images.
void downsample_internal_with_mask(const float* image_in, const float * mask_in,
		float* image_out, float *mask_out, uint size_x, uint size_y, 
		uint sample_size_x, uint sample_size_y);

void upsample_internal(const float* in, float* out, uint size_x, uint size_y, uint sample_size_x, uint sample_size_y);

int get_sample_size(uint size, float factor);
void separate_convolution(const float *in, float *out, int size_x, int size_y, const float *filter_x, const float *filter_y, int filter_x_size, int filter_y_size);

float bilinear_interpolation(const float *input, int size_x, int size_y, float x, float y);
inline float nearest_interpolation(const float *input, int size_x, int size_y, float x, float y);

inline int neumann_boundary_condition(int x, int size);
inline int periodic_boundary_condition(int x, int size);
inline int symmetric_boundary_condition_A(int x, int size);
inline int symmetric_boundary_condition_B(int x, int size);
inline int cutting_boundary_condition(int x, int size, bool &is_out);

} // namespace _Details

} // namespace Sampling


#endif /* SAMPLING_H_ */
