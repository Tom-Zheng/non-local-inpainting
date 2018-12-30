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

#include "sampling.h"
#include <cmath>
#include <cstdlib>

namespace Sampling {

Image<float> downsample(FixedImage<float> in, float factor)
{
	// check if factor is valid
	if (factor <= 0 || factor >= 1.0) {
		return Image<float>();
	}

	// size and channels of downsampled image
	int sample_size_x = _Details::get_sample_size(in.get_size_x(), factor);
	int sample_size_y = _Details::get_sample_size(in.get_size_y(), factor);
	uint number_of_channels = in.get_number_of_channels();

	// allocate output
	Image<float> out(sample_size_x, sample_size_y, number_of_channels);

	// perform downsampling by gaussian filtering and bilinear interp
	if (number_of_channels == 1) {
		_Details::downsample_internal(in.raw(), out.raw(), in.get_size_x(), in.get_size_y(), sample_size_x, sample_size_y);
	} else {
		// process each channel separatedly
		const float *in_data = in.raw();
		float *out_data = out.raw();

		float *channel_in = new float[in.get_size_x() * in.get_size_y()];
		float *channel_out = new float[sample_size_x * sample_size_y];
		for (uint ch = 0; ch < number_of_channels; ch++) {
			// extract single channel from the input data
			for (uint y = 0; y < in.get_size_y(); y++) {
				for (uint x = 0; x < in.get_size_x(); x++) {
					int index = y * in.get_size_x() + x;
					channel_in[index] = in_data[index * number_of_channels + ch];
				}
			}

			// downsample single channel
			_Details::downsample_internal(channel_in, channel_out, in.get_size_x(), in.get_size_y(), sample_size_x, sample_size_y);

			// copy processed channel into the output data
			for (int y = 0; y < sample_size_y; y++) {
				for (int x = 0; x < sample_size_x; x++) {
					int index = y * sample_size_x + x;
					out_data[index * number_of_channels + ch] = channel_out[index];
				}
			}
		}
		delete[] channel_in;
		delete[] channel_out;
	}

	return out;
}


Mask downsample(FixedMask in, float factor, float threshold)
{
	// check if factor is valid
	if (factor <= 0 || factor >= 1.0) {
		return Mask();
	}

	// convert mask from binary to float
	float *float_mask = new float[in.get_size_x() * in.get_size_y()]();
	FixedMask::iterator it;
	for (it = in.begin(); it != in.end(); ++it) {
		int index = it->y * in.get_size_x() +  it->x;
		float_mask[index] = 1.0;
	}

	// size and channels of downsampled image
	int sample_size_x = _Details::get_sample_size(in.get_size_x(), factor);
	int sample_size_y = _Details::get_sample_size(in.get_size_y(), factor);
	float *sampled_float_mask = new float[sample_size_x * sample_size_y]();

	// perform downsampling by gaussian filtering and bilinear interp
	_Details::downsample_internal(float_mask, sampled_float_mask, in.get_size_x(), in.get_size_y(), sample_size_x, sample_size_y);

	// float to binary by thresholding with specified threshold
	threshold = fmax(0.0, fmin(1.0, threshold));
	Mask out(sample_size_x, sample_size_y);
	float *p_sampled_float_mask = sampled_float_mask;
	for (int y = 0; y < sample_size_y; y++) {
		for (int x = 0; x < sample_size_x; x++, p_sampled_float_mask++) {
			if (*p_sampled_float_mask > threshold) {
				out.mask(x, y);
			}
		}
	}

	delete[] sampled_float_mask;
	delete[] float_mask;

	return out;
}


std::pair<Image<float>, Mask> downsample_with_mask(FixedImage<float> image_in,
		Mask mask_in, float factor, float threshold, bool mask_val)
{
	typedef std::pair<Image<float>, Mask> OutputPair;

	// check if factor is valid
	if (factor <= 0 || factor >= 1.0) {
		return OutputPair(Image<float>(), Mask());
	}

	// check if image and mask have same size
	if (image_in.get_size() != mask_in.get_size()) {
		return OutputPair(Image<float>(), Mask());
	}

	// convert mask from binary to float
	float *float_mask = new float[mask_in.get_size_x() * mask_in.get_size_y()]();
	float *p_float_mask = float_mask;
	for (uint y = 0; y < mask_in.get_size_y(); y++) {
		for (uint x = 0; x < mask_in.get_size_x(); x++, p_float_mask++) {
			*p_float_mask = (mask_val == mask_in(x,y)) ? 1.0 : 0.0;
		}
	}

	// size and channels of downsampled image
	int sample_size_x = _Details::get_sample_size(image_in.get_size_x(), factor);
	int sample_size_y = _Details::get_sample_size(image_in.get_size_y(), factor);
	uint number_of_channels = image_in.get_number_of_channels();

	// allocate outputs and intermediate mask buffer
	OutputPair outputs(Image<float>(sample_size_x, sample_size_y, number_of_channels),
	                   Mask(sample_size_x, sample_size_y));
	Image<float> image_out = outputs.first;
	float *sampled_float_mask = new float[sample_size_x * sample_size_y]();
	Mask mask_out = outputs.second;

	// perform downsampling by gaussian filtering and bilinear interp
	if (number_of_channels == 1) {
		_Details::downsample_internal_with_mask(image_in.raw(), float_mask,
		                      image_out.raw(), sampled_float_mask, 
		                      image_in.get_size_x(), image_in.get_size_y(),
		                      sample_size_x, sample_size_y);
	} else {
		// process each channel separatedly
		const float *in_data = image_in.raw();
		float *out_data = image_out.raw();

		float *channel_in = new float[image_in.get_size_x() * image_in.get_size_y()];
		float *channel_out = new float[sample_size_x * sample_size_y];
		for (uint ch = 0; ch < number_of_channels; ch++) {
			// extract single channel from the input data
			for (uint y = 0; y < image_in.get_size_y(); y++) {
				for (uint x = 0; x < image_in.get_size_x(); x++) {
					int index = y * image_in.get_size_x() + x;
					channel_in[index] = in_data[index * number_of_channels + ch];
				}
			}

			// downsample single channel
			_Details::downsample_internal_with_mask(channel_in, float_mask,
			                     channel_out, sampled_float_mask, 
			                     image_in.get_size_x(), image_in.get_size_y(),
			                     sample_size_x, sample_size_y);

			// copy processed channel into the output data
			for (int y = 0; y < sample_size_y; y++) {
				for (int x = 0; x < sample_size_x; x++) {
					int index = y * sample_size_x + x;
					out_data[index * number_of_channels + ch] = channel_out[index];
				}
			}
		}

		delete[] channel_in;
		delete[] channel_out;
	}

	// float to binary by thresholding with specified threshold
	threshold = fmax(0.0, fmin(1.0, threshold));
	float *p_sampled_float_mask = sampled_float_mask;
	for (int y = 0; y < sample_size_y; y++) {
		for (int x = 0; x < sample_size_x; x++, p_sampled_float_mask++) {
			mask_out(x,y) = mask_val ? (*p_sampled_float_mask > threshold) : 
			                           (*p_sampled_float_mask < threshold) ; 
		}
	}

	return outputs;
}


Image<float> upsample(FixedImage<float> in, Shape size)
{
	int sample_size_x = size.size_x;
	int sample_size_y = size.size_y;
	uint number_of_channels = in.get_number_of_channels();

	Image<float> out(sample_size_x, sample_size_y, number_of_channels);

	if (number_of_channels == 1) {
		_Details::upsample_internal(in.raw(), out.raw(), in.get_size_x(), in.get_size_y(), sample_size_x, sample_size_y);
	} else {
		const float *in_data = in.raw();
		float *out_data = out.raw();

		float *channel_in = new float[in.get_size_x() * in.get_size_y()];
		float *channel_out = new float[sample_size_x * sample_size_y];
		for (uint ch = 0; ch < number_of_channels; ch++) {
			// extract single channel from the input data
			for (uint y = 0; y < in.get_size_y(); y++) {
				for (uint x = 0; x < in.get_size_x(); x++) {
					int index = y * in.get_size_x() + x;
					channel_in[index] = in_data[index * number_of_channels + ch];
				}
			}

			// upsample single channel
			_Details::upsample_internal(channel_in, channel_out, in.get_size_x(), in.get_size_y(), sample_size_x, sample_size_y);

			// copy processed channel into the output data
			for (int y = 0; y < sample_size_y; y++) {
				for (int x = 0; x < sample_size_x; x++) {
					int index = y * sample_size_x + x;
					out_data[index * number_of_channels + ch] = channel_out[index];
				}
			}
		}
		delete[] channel_in;
		delete[] channel_out;
	}

	return out;
}

namespace _Details {

void downsample_internal(const float* in, float* out, uint size_x, uint size_y, uint sample_size_x, uint sample_size_y)
{
	// downsampling factor
	float factor_x = (float)size_x / sample_size_x;
	float factor_y = (float)size_y / sample_size_y;

	// adjust smoothing filter
	float sigma_x = 0.62 * sqrt( pow(factor_x, 2) - 1 );
	float sigma_y = 0.62 * sqrt( pow(factor_y, 2) - 1 );
	//int kernel_size_x = 2 * round(1.5 * sigma_x) + 1;
	//int kernel_size_y = 2 * round(1.5 * sigma_y) + 1;
	int kernel_size_x = 9; // HARDCODED to speed up pyramid
	int kernel_size_y = 9;

	// horizontal and vertical 1D Gaussian kernels
	Image<float> gaussian_x = GaussianWeights::calculate_1d(kernel_size_x, sigma_x);
	const float *filter_x = gaussian_x.raw();
	Image<float> gaussian_y = GaussianWeights::calculate_1d(kernel_size_y, sigma_y);
	const float *filter_y = gaussian_y.raw();

	// smooth image using separable convolution - result will be stored in `buffer`
	float *buffer = new float[size_x * size_y];
	_Details::separate_convolution(in, buffer, size_x, size_y,
	                               filter_x, filter_y, kernel_size_x, kernel_size_y);

	// subsample smoothed image using bilinear interp

	// we center the coarse grid over the fine grid by adding this offset the
	// scaled coarse grid coordinates:
	float min_x = factor_x / 2.0 - 0.5;
	float min_y = factor_y / 2.0 - 0.5;

	/* NOTE: The sampling grid
	 *
	 * x coords in the fine   grid are (0, 1, ...,        size_x - 1)
	 * x coords in the coarse grid are (0, 1, ..., sample_size_x - 1) * factor_x + min_x
	 *
	 * min_x is chosen such that the coarse grid is centered over the
	 * fine grid: the center of both grids coincide. */

	for(uint i = 0; i < sample_size_y; i++) {
		float y = min_y + i * factor_y;
		for(uint j = 0; j < sample_size_x; j++) {
			float x = min_x + j * factor_x;

			float value = _Details::bilinear_interpolation(buffer, size_x, size_y, x, y);
			out[i * sample_size_x + j] = value;
		}
	}

	// free memory
	delete [] buffer;
}

void downsample_internal_with_mask(const float* img_in, const float * msk_in,
		float* img_out, float *msk_out, uint size_x, uint size_y, 
		uint sample_size_x, uint sample_size_y)
{
	// downsampling factor
	float factor_x = (float)size_x / sample_size_x;
	float factor_y = (float)size_y / sample_size_y;

	// multiply image by the mask
	float *img_x_msk = new float[size_x * size_y];
	float *p_img_x_msk = img_x_msk;
	const float *p_img = img_in, *p_msk = msk_in;
	for (uint y = 0; y < size_y; y++) {
		for (uint x = 0; x < size_x; x++, p_img_x_msk++, p_img++, p_msk++) {
			*p_img_x_msk = *p_msk ? *p_img : 0.;
		}
	}

	// filter msk and masked image
	float *msk_filtered = new float[size_x * size_y];
	float *img_x_msk_filtered = new float[size_x * size_y];
	{
		// adjust smoothing filter
		float sigma_x = 0.62 * sqrt( pow(factor_x, 2) - 1 );
		float sigma_y = 0.62 * sqrt( pow(factor_y, 2) - 1 );
		//int kernel_size_x = 2 * round(1.5 * sigma_x) + 1;
		//int kernel_size_y = 2 * round(1.5 * sigma_y) + 1;
		int kernel_size_x = 9; // HARDCODED to speed up pyramid
		int kernel_size_y = 9;

		// horizontal and vertical 1D Gaussian kernels
		Image<float> gaussian_x = GaussianWeights::calculate_1d(kernel_size_x, sigma_x);
		Image<float> gaussian_y = GaussianWeights::calculate_1d(kernel_size_y, sigma_y);
		const float *filter_x = gaussian_x.raw();
		const float *filter_y = gaussian_y.raw();

		// filter mask using separable convolution
		_Details::separate_convolution(msk_in, msk_filtered, size_x, size_y,
				filter_x, filter_y, kernel_size_x, kernel_size_y);

		// filter masked image using separable convolution
		_Details::separate_convolution(img_x_msk, img_x_msk_filtered, size_x, size_y,
				filter_x, filter_y, kernel_size_x, kernel_size_y);
	}

	// subsample smoothed images using bilinear interp
	{
		// we center the coarse grid over the fine grid by adding this offset the
		// scaled coarse grid coordinates:
		float min_x = factor_x / 2.0 - 0.5;
		float min_y = factor_y / 2.0 - 0.5;

		/* NOTE: The sampling grid
		 *
		 * x coords in the fine   grid are (0, 1, ...,        size_x - 1)
		 * x coords in the coarse grid are (0, 1, ..., sample_size_x - 1) * factor_x + min_x
		 *
		 * min_x is chosen such that the coarse grid is centered over the
		 * fine grid: the center of both grids coincide. */

		using _Details::bilinear_interpolation;
		float *p_img_out = img_out;
		float *p_msk_out = msk_out;
		for(uint i = 0; i < sample_size_y; i++) {
			float y = min_y + i * factor_y;
			for(uint j = 0; j < sample_size_x; j++, p_img_out++, p_msk_out++) {
				float x = min_x + j * factor_x;

				float msk_value = bilinear_interpolation(msk_filtered      , size_x, size_y, x, y);
				float img_value = bilinear_interpolation(img_x_msk_filtered, size_x, size_y, x, y);

				*p_img_out = (msk_value > 1e-6) ? img_value / msk_value : 0.f ;
				*p_msk_out = msk_value;
			}
		}
	}

	// free memory
	delete [] img_x_msk;
	delete [] img_x_msk_filtered;
	delete [] msk_filtered;
}

void upsample_internal(const float* in, float* out, uint size_x, uint size_y, uint sample_size_x, uint sample_size_y)
{
	float factor_x = (float) sample_size_x / size_x;
	float factor_y = (float) sample_size_y / size_y;

	// bounding box
	float min_x = factor_x / 2.0 - 0.5;
	float min_y = factor_y / 2.0 - 0.5;

	// set samples in output image
	for(uint i = 0; i < sample_size_y; i++) {
		float y = (i - min_y) / factor_y;
		for(uint j = 0; j < sample_size_x; j++) {
			float x = (j - min_x) / factor_x;
			float value = _Details::nearest_interpolation(in, size_x, size_y, x, y);
			out[i * sample_size_x + j] = value;
		}
	}
}


int get_sample_size(uint size, float factor)
{
	return (int)((float) size * factor + 0.5);
}


void separate_convolution(const float *in, float *out, int size_x, int size_y, const float *filter_x, const float *filter_y, int filter_x_size, int filter_y_size)
{
	// initialize temporal buffer
	float *buffer;
	buffer = new float[size_x * size_y];

	float sum;
	int id;

	// convolution along x axis
	int radius = (filter_x_size - 1) / 2;

	for (int y = 0; y < size_y; y++) {
		for (int x = 0;x < size_x; x++) {
			sum = 0.0;

			for (int i = filter_x_size - 1; i >= 0; i--) {
				id = x + radius - i;
				id = symmetric_boundary_condition_A(id, size_x);

				sum += filter_x[i] * in[y * size_x + id];
			}

			buffer[y * size_x + x] = sum;
		}
	}

	// convolution along y axis
	radius = (filter_y_size - 1) / 2;

	for (int y = 0;y < size_y; y++) {
		for (int x = 0; x < size_x; x++) {
			sum = 0.0;

			for (int i = filter_y_size - 1; i >= 0; i--) {
				id = y + radius - i;
				id = symmetric_boundary_condition_A(id, size_y);

				sum += filter_y[i] * buffer[id * size_x + x];
			}

			out[y * size_x + x] = sum;
		}
	}

	// free memory
	delete [] buffer;
}


float bilinear_interpolation(const float *input, int size_x, int size_y, float x, float y)
{
	int id_x = floor(x);
	int id_y = floor(y);

	int ids_x[2], ids_y[2];

	// apply the appropriate boundary conditions
	ids_x[0] = neumann_boundary_condition(id_x, size_x);
	ids_y[0] = neumann_boundary_condition(id_y, size_y);
	ids_x[1] = neumann_boundary_condition(id_x + 1, size_x);
	ids_y[1] = neumann_boundary_condition(id_y + 1, size_y);

	const float p11 = input[ ids_x[0] + size_x * ids_y[0] ];

	float a = x - id_x;
	float b = y - id_y;
	float result = 0.0;

	// interpolate
	if ((!a) || (!b)) {
		if ((!a) && (!b)) {
			result = p11;
		} else {
			if (!a) {
				const float p12 = input[ ids_x[0] + size_x * ids_y[1] ];
				result = (1.0 - b) * p11 + b * p12;
			} else {
				const float p21 = input[ ids_x[1] + size_x * ids_y[0] ];
				result = (1.0 - a) * p11 + a * p21;
			}
		}
	} else {
		float a_1 = 1.0 - a;
		float b_1 = 1.0 - b;
		const float p21 = input[ ids_x[1] + size_x * ids_y[0] ];
		const float p12 = input[ ids_x[0] + size_x * ids_y[1] ];
		const float p22 = input[ ids_x[1] + size_x * ids_y[1] ];
		result = b_1 * (a_1 * p11 + a * p21) + b * (a_1 * p12 + a * p22);
	}

	return result;
}


inline float nearest_interpolation(const float *input, int size_x, int size_y, float x, float y)
{
	int id_x = floor(x + 0.5);	// NOTE: this is 'round'
	int id_y = floor(y + 0.5);	// NOTE: this is 'round'

	// apply the appropriate boundary conditions
	id_x = neumann_boundary_condition(id_x, size_x);
	id_y = neumann_boundary_condition(id_y, size_y);

	float result = input[id_x + size_x * id_y];

	return result;
}

/**
 * For size = 4 : | 0 0 0 0 | 0 1 2 3 | 3 3 3 3 |
 */
inline int neumann_boundary_condition(int x, int size)
{
	return max(0, min(size - 1, x));
}


/**
 * For size = 4 : | 0 1 2 3 | 0 1 2 3 | 0 1 2 3 |
 */
inline int periodic_boundary_condition(int x, int size)
{
	if(x < 0) {
		x = (size - std::abs(x) % size) % size;
	} else if(x >= size) {
		x = x % size;
	}

	return x;
}


/**
 * For size = 4 : | 3 2 1 0 | 0 1 2 3 | 3 2 1 0 |
 */
inline int symmetric_boundary_condition_A(int x, int size)
{
	while ((x < 0) || (x >= size)) {
		if (x < 0) {
			x = -x - 1;
		}
		if (x >= size) {
			x = 2 * size - x -1;
		}
	}

	return x;
}


/**
 * For size = 4 : | 2 3 2 1 | 0 1 2 3 | 2 1 0 1 |
 */
inline int symmetric_boundary_condition_B(int x, int size)
{
	if(x < 0 || x >= size) {
		const int border = size - 1;
		const int abs_x = std::abs(x);

		x = ( (abs_x / border) % 2 ) ? border - ( abs_x % border ) : abs_x % border;
	}

	return x;
}

/**
 * For size = 4 : | -4 -3 -2 -1 | 0 1 2 3 | 4 5 6 7 |
 */
inline int cutting_boundary_condition(int x, int size, bool &is_out)
{
	if(x < 0 || x >= size) {
		is_out = true;
	}

	return x;
}

} // namespace _Detailes

} // namespace Sampling
