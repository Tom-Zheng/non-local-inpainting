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

#include "image_inpainting.h"
#include <stdexcept>

extern "C" {
   // defined in 3rdparty/simpois/simpois.c
   void poisson_solver_separable_simplest(float *out, float *in, float *dat,
        int w, int h, int pd, int niter);
}

const float ImageInpainting::MASK_SAMPLING_THRESHOLD = 0.4f;

ImageInpainting::ImageInpainting()
{
	_patch_match = 0;
	_image_updating = 0;
	_iterations_amount = 5;
	_tolerance = 0;
	_scales_amount = 3;
	_subsampling_rate = 0.5;
	_confidence_decay_time = 5.0;
	_confidence_asymptotic_value = 0.1;
	_initialization_type = InitPoisson;

	_keep_intermediate = false;
}


ImageInpainting::ImageInpainting(int iterations_amount,
								 float tolerance,
								 int scales_amount,
								 float subsampling_rate,
								 float confidence_decay_time,
								 float confidence_asymptotic_value,
								 InitType initialization_type)
{
	_patch_match = 0;
	_image_updating = 0;
	_iterations_amount = iterations_amount;
	_tolerance = tolerance;
	_scales_amount = scales_amount;
	_subsampling_rate = subsampling_rate;
	_confidence_decay_time = confidence_decay_time;
	_confidence_asymptotic_value = confidence_asymptotic_value;
	_initialization_type = initialization_type;

	_keep_intermediate = false;
}


/**
 * Computes the multiscale image inpainting
 *
 * @param in Image to be inpainted
 * @param mask Inpainting domain
 */
Image<float> ImageInpainting::process(FixedImage<float> in, FixedMask mask)
{
	_nnf_pyramid.clear();
	_original_image_pyramid.clear();
	_image_pyramid.resize(_scales_amount);
	_mask_pyramid.resize(_scales_amount);
	_image_pyramid[0] = Image<float>(in);	// this will be changed, thus we make a copy
	_mask_pyramid[0] = mask;

#ifdef DBG_OUTPUT
	cm_ind = 0;
#endif

	printf("\tbuilding image pyramid\n");
	for (int i = 1; i < _scales_amount; i++) {
		if (_initialization_type != InitNone) {
			std::pair<Image<float>, FixedMask> downscaled_pair;
			downscaled_pair = Sampling::downsample_with_mask(
					_image_pyramid[i - 1], _mask_pyramid[i - 1],
					_subsampling_rate, MASK_SAMPLING_THRESHOLD, false);

			_image_pyramid[i] = downscaled_pair.first;
			_mask_pyramid[i] = downscaled_pair.second;
		} else {
			_image_pyramid[i] = Sampling::downsample(_image_pyramid[i - 1], _subsampling_rate);
			_mask_pyramid[i] = Sampling::downsample(_mask_pyramid[i - 1], _subsampling_rate, MASK_SAMPLING_THRESHOLD);
		}

#ifdef DBG_OUTPUT
		IOUtility::write_rgb_image (IOUtility::compose_file_name("dbg_level", _scales_amount - i - 1, "png"), IOUtility::lab_to_rgb(_image_pyramid[i]));
		IOUtility::write_mono_image(IOUtility::compose_file_name("dbg_level_mask", _scales_amount - i - 1, "pgm"), IOUtility::mask_to_greyscale(_mask_pyramid[i]));
#endif
	}

	// keep original image pyramid, if needed
	if (_keep_intermediate) {
		vector<Image<float> >::iterator it;
		for (it = _image_pyramid.begin(); it != _image_pyramid.end(); ++it) {
			_original_image_pyramid.push_back(it->clone());
		}
	}
	
	// initialize at the coarsest level, if needed
	if (_initialization_type != InitNone) {
        switch (_initialization_type) {
            case InitAvg: {
		        float *initial_color = get_initial_color(_image_pyramid.back(), _mask_pyramid.back(), true);
		        initialize_with_color(_image_pyramid.back(), _mask_pyramid.back(), initial_color);
		        delete initial_color;
                break;
            }
            case InitBlack: { // get_initial_color with parameter false returns 0's
		        float *initial_color = get_initial_color(_image_pyramid.back(), _mask_pyramid.back(), false);
		        initialize_with_color(_image_pyramid.back(), _mask_pyramid.back(), initial_color);
		        delete initial_color;
                break;
            }
            case InitPoisson: { 
                // prepare data for interfacing external "simple poisson solver"
                int w  = _image_pyramid.back().get_size_x(); 
                int h  = _image_pyramid.back().get_size_y(); 
                int ch = _image_pyramid.back().get_number_of_channels();
                int sz = w*h*ch;
                float * in   = new float[sz];
                float * out  = new float[sz];

                for (int c = 0; c < ch; c++) 
                    for (int j = 0; j < h; j++)
                        for (int i = 0; i < w; i++) {
                            if (_mask_pyramid.back().test(i,j)) 
                                in[i + j*w + c*w*h] = NAN;  // solver uses NAN to signal mask
                            else
                                in[i + j*w + c*w*h] = _image_pyramid.back()(i,j,c);
                        }

                poisson_solver_separable_simplest(out, in, NULL, w, h, ch, 10);

                for (int c = 0; c < ch; c++)  // copy the result back to image_pyramid
                    for (int j = 0; j < h; j++)
                        for (int i = 0; i < w; i++) 
                            _image_pyramid.back()(i,j,c) = out[i + j*w + c*w*h];

                delete[] out;
                delete[] in;
                break;
            }
        }
	}

#ifdef DBG_OUTPUT
	IOUtility::write_rgb_image(IOUtility::compose_file_name("dbg_coarsest_initialized.png"), IOUtility::lab_to_rgb(_image_pyramid.back()));
#endif

	Image<float> confidence_mask = calculate_confidence_mask(_mask_pyramid.back(), _confidence_decay_time, _confidence_asymptotic_value);

	Image<Point> nnf;

	printf("\tinpainting scale %d\n", _scales_amount);
	inpaint_internal(_image_pyramid.back(), _mask_pyramid.back(), confidence_mask, nnf, _tolerance);

#ifdef DBG_OUTPUT
	IOUtility::write_rgb_image(IOUtility::compose_file_name("dbg_inpainted_0.png"), IOUtility::lab_to_rgb(_image_pyramid.back()));
	IOUtility::write_mono_image(IOUtility::compose_file_name("dbg_confidence_mask_0.png"), IOUtility::probability_to_greyscale(confidence_mask));
#endif

	for (int i = _scales_amount - 2; i >= 0; i--) {
		Image<float> confidence_mask = calculate_confidence_mask(_mask_pyramid[i], _confidence_decay_time, _confidence_asymptotic_value);

		nnf = propagate_weights(_image_pyramid[i], _mask_pyramid[i], confidence_mask, _image_pyramid[i + 1], _mask_pyramid[i + 1]);

#ifdef DBG_OUTPUT
		IOUtility::write_rgb_image(IOUtility::compose_file_name("dbg_propagated", _scales_amount - i - 1, "png"), IOUtility::lab_to_rgb(_image_pyramid[i]));
		IOUtility::write_mono_image(IOUtility::compose_file_name("dbg_confidence_mask", _scales_amount - i - 1, "png"), IOUtility::probability_to_greyscale(confidence_mask));
#endif
		printf("\tinpainting scale %d\n", i+1);
		inpaint_internal(_image_pyramid[i], _mask_pyramid[i], confidence_mask, nnf, _tolerance);
#ifdef DBG_OUTPUT
		IOUtility::write_rgb_image(IOUtility::compose_file_name("dbg_inpainted", _scales_amount - i - 1, "png"), IOUtility::lab_to_rgb(_image_pyramid[i]));
#endif
	}
	
	return _image_pyramid.front();
}


/**
 * @param size_ratio Size of the image at the coarsest scale to its size at the finest scale (coarsest/finest)
 * @param scales_amount Number of levels in the pyramid
 */
float ImageInpainting::calculate_subsampling_rate(float size_ratio, unsigned int scales_amount)
{
	if (size_ratio < 0.0 || size_ratio > 1.0) {
		return 0.0;
	}

	return exp(log(size_ratio) / (scales_amount - 1));
}


int ImageInpainting::get_iterations_amount()
{
	return _iterations_amount;
}


void ImageInpainting::set_iterations_amount(int amount)
{
	_iterations_amount = amount;
}


float ImageInpainting::get_tolerance()
{
	return _tolerance;
}


void ImageInpainting::set_tolerance(float value)
{
	_tolerance = value;
}


int ImageInpainting::get_scales_amount()
{
	return _scales_amount;
}


void ImageInpainting::set_scales_amount(int amount)
{
	_scales_amount = amount;
}


float ImageInpainting::get_subsampling_rate()
{
	return _subsampling_rate;
}


void ImageInpainting::set_subsampling_rate(float rate)
{
	_subsampling_rate = rate;
}


float ImageInpainting::get_confidence_decay_time()
{
	return _confidence_decay_time;
}


void ImageInpainting::set_confidence_decay_time(float value)
{
	_confidence_decay_time = value;
}


float ImageInpainting::get_confidence_asymptotic_value()
{
	return _confidence_asymptotic_value;
}


void ImageInpainting::set_confidence_asymptotic_value(float value)
{
	_confidence_asymptotic_value = value;
}


PatchMatch* ImageInpainting::get_weights_updating()
{
	return _patch_match;
}


void ImageInpainting::set_weights_updating(PatchMatch *patch_match)
{
	_patch_match = patch_match;
}


AImageUpdating* ImageInpainting::get_image_updating()
{
	return _image_updating;
}


void ImageInpainting::set_image_updating(AImageUpdating *image_updating)
{
	_image_updating = image_updating;
}


/**
 * Specifies whether intermediate data ('input image pyramid' and 'nnf pyramid')
 * should be stored, or not. Note: 'output image pyramid' and 'mask pyramid' will
 * be stored anyway, since they are used in computations.
 *
 * @param value True to keep intermediate data.
 */
void ImageInpainting::keep_intermediate(bool value)
{
	_keep_intermediate = value;
}


vector<Image<float> > ImageInpainting::get_input_image_pyramid()
{
	return _original_image_pyramid;
}


vector<Image<float> > ImageInpainting::get_output_image_pyramid()
{
	return _image_pyramid;
}


vector<Mask> ImageInpainting::get_mask_pyramid()
{
	return _mask_pyramid;
}


vector<Image<Point> > ImageInpainting::get_nnf_pyramid()
{
	return _nnf_pyramid;
}

/* Private */

/**
 * Updates the image by computing the single scale inpainting
 *
 * @param image Image to be inpainted
 * @param mask Binary mask of the inpainting domain
 * @param confidence_mask Confidence values from a range [0.0, 1.0] for every point on the image
 * @param tolerance Stopping criteria
 */
void ImageInpainting::inpaint_internal(Image<float> image,
									   FixedMask inpainting_domain,
									   FixedImage<float> confidence_mask,
									   FixedImage<Point> initial_nnf,
									   float tolerance)
{
	Shape patch_size = _image_updating->get_patch_size();

	Image<float> original_image = image.clone();

	Mask extended_inpainting_domain = get_extended_domain(inpainting_domain, patch_size);

	Mask target_mask = extended_inpainting_domain.clone();			// make an explicit copy
	Mask source_mask = extended_inpainting_domain.clone_invert();	// make an explicit copy and invert

	// modify the source_mask to consider the incomplete forward gradients that require
	// an extra unmasked pixel to be computed (definition of O^c_e section 3.2 of the IPOL article)
	// A less deliacte alternative is to erode the whole source_mask
	// Mask source_mask = get_extended_domain(target_mask, Shape(3,3)).clone_invert();
	bool is_source_mask_empty = true;
	for (uint j=0; j< source_mask.get_size_y(); j++) {
		for (uint i=0; i< source_mask.get_size_x(); i++) {
			if (source_mask.test(i,j)==true) {
				if      (source_mask.test(i+1, j)==false) source_mask.unmask(i, j);
				else if (source_mask.test(i, j+1)==false) source_mask.unmask(i, j);
				else is_source_mask_empty = false;
			}
		}
	}
	// add margin at the border of source_mask and target_mask
	add_margin(source_mask, target_mask, patch_size.size_x / 2);

	// check if source mask is empty
	if (is_source_mask_empty) {
		throw std::runtime_error("ERROR: Empty source mask (no complete patches to copy from. This may happen due to a too big inpainting domain, too big patch, or too much downscaling)");
	}

	Image<Point> nnf = initial_nnf;
	double total_difference = numeric_limits<double>::max();
	int i = 0;
	for (i = 0; i < _iterations_amount && total_difference > tolerance; i++) {
		// update weights (find nearest neighbours field)
		nnf = _patch_match->calculate(image, source_mask, image, target_mask, nnf);

		// update image
		total_difference = _image_updating->update(image, original_image, inpainting_domain, extended_inpainting_domain, nnf, confidence_mask);
		// DEBUG
		if (i % 10 == 0)
			printf("\t\t\tIter:%d, Err: %f\n", i, total_difference);
#ifdef DBG_OUTPUT
		if (i % 10 == 0)
			IOUtility::write_rgb_image(IOUtility::compose_file_name("dbg_inpainted", cm_ind, i, "png"), IOUtility::lab_to_rgb(image));
#endif
	}
	printf("\t\tinpainting ended after %d iterations\n", i);

	// keep nnf, if needed
	if (_keep_intermediate) {
		_nnf_pyramid.push_back(nnf);
	}

#ifdef DBG_OUTPUT
	cm_ind++;
#endif
}


/**
 * Updates the image by filling all the pixels given by the mask with the given color value
 *
 * @param image Image to be updated
 * @param mask Binary mask which specifies pixels to be filled
 * @param color Color to assign to masked pixels
 */
void ImageInpainting::initialize_with_color(Image<float> image, FixedMask mask, float *color)
{
	FixedMask::iterator it;
	for (it = mask.begin(); it != mask.end(); ++it) {
		for (uint channel = 0; channel < image.get_number_of_channels(); channel++) {
			image(*it, channel) = color[channel];
		}
	}
}


/**
 * Calculates initial color for initialization
 * @note we use a dedicated method to get initial color to easily experiment with different colors
 *
 * @param image Image to be initialized
 * @param mask Binary mask which specifies the inpainting domain
 * @param average Flag which specifies that the average color outside the inpainting domain should be calculated. Otherwise this function returns black.
 */
float* ImageInpainting::get_initial_color(Image<float> image,
										  FixedMask mask,
										  bool average)
{
	uint channels = image.get_number_of_channels();
	float* color = new float[channels]();

	// return black
	if (!average) {
		for (uint ch = 0; ch < channels; ch++) {
			color[ch] = 0;
		}
		return color;
	}

	// return average of known region
	for (uint ch = 0; ch < channels; ch++) {
		float sum = 0.0f;
		int count = 0;
		for (uint y = 0; y < image.get_size_y(); y++) {
			for (uint x = 0; x < image.get_size_x(); x++) {
				if (!mask(x, y)) {
					sum += image(x, y, ch);
					count++;
				}
			}
		}

		if (count > 0) {
			color[ch] = sum / (float)count;
		}
	}

	return color;
}


/**
 * Updates the image by propagating information from the lower level of the pyramid into the inpainting domain.
 *
 * @param upper_level Image to be updated (current level of image pyramid)
 * @param upper_inpainting_domain Corresponding inpainting domain (current level of mask pyramid)
 * @param upper_confidence_mask Corresponding confidence values
 * @param lower_level Image to take information from (lower level of image pyramid)
 * @param lower_inpainting_domain Corresponding inpainting domain (lower level of mask pyramid)
 */
Image<Point> ImageInpainting::propagate_weights(Image<float> upper_level,
										FixedMask upper_inpainting_domain,
										FixedImage<float> upper_confidence_mask,
										FixedImage<float> lower_level,
										FixedMask lower_inpainting_domain)
{
	Shape patch_size = _image_updating->get_patch_size();

	// prepare the masks
	Mask lower_target_mask = get_extended_domain(lower_inpainting_domain, patch_size);
	Mask lower_source_mask = lower_target_mask.clone_invert();

	// add margin at the border of lower_target_mask and lower_source_mask
	int half_patch_side = _image_updating->get_patch_size().size_x / 2;
	add_margin(lower_source_mask, lower_target_mask, half_patch_side);

	// calculate NNF using image and inpainting domain mask from the lower level
	Image<Point> nnf = _patch_match->calculate(lower_level, lower_source_mask, lower_level, lower_target_mask);

	// scale NNF
	float scale_x = (float)upper_level.get_size_x() / lower_level.get_size_x();
	float scale_y = (float)upper_level.get_size_y() / lower_level.get_size_y();

	// upsample scaled NNF using nearest-neighbor interpolation
	Image<Point> scaled_nnf(upper_level.get_size());
	scaled_nnf.fill(Point(-1, -1));

	FixedMask::iterator it;
	for (it = lower_inpainting_domain.begin(); it != lower_inpainting_domain.end(); ++it) {
		Point neighbor = nnf(it->x, it->y);
		if (neighbor.x >= 0) {
			neighbor.x = round((float)neighbor.x * scale_x);
			neighbor.y = round((float)neighbor.y * scale_y);
			int x = round((float)it->x * scale_x);
			int y = round((float)it->y * scale_y);

			scaled_nnf(x, y) = neighbor;
		}
	}

	// upsample the lower level image
	Image<float> upsampled_image = Sampling::upsample(lower_level, upper_level.get_size());

	// copy upsampled information inside the inpainting domain to initialize the upper level
	for (it = upper_inpainting_domain.begin(); it != upper_inpainting_domain.end(); ++it) {
		for (uint channel = 0; channel < upsampled_image.get_number_of_channels(); channel++) {
		 	 upper_level(*it, channel) = upsampled_image(*it, channel);
		}
	}

	Mask extended_upper_inpainting_domain = get_extended_domain(upper_inpainting_domain, patch_size);

	// refine NNF using the scaled_nnf and the upsampled image
	Mask upper_target_mask = extended_upper_inpainting_domain.clone();
	Mask upper_source_mask = extended_upper_inpainting_domain.clone_invert();

	// add margin at the border of lower_target_mask and lower_source_mask
	add_margin(upper_source_mask, upper_target_mask, half_patch_side);

	Image<Point> refined_nnf = _patch_match->calculate(upsampled_image, upper_source_mask, upsampled_image, upper_target_mask, scaled_nnf);

	// use scaled NNF and rough initialization to propagate information to the upper level
	_image_updating->update(upper_level, upper_level, upper_inpainting_domain, extended_upper_inpainting_domain, refined_nnf, upper_confidence_mask);

	return refined_nnf;
}


/**
 * Computes the mask containing centers of all patches, intersecting the given inpainting domain
 *
 * @param inpainting_domain Inpainting domain to be extended
 * @param patch_size Size of patches
 */
Mask ImageInpainting::get_extended_domain(FixedMask inpainting_domain, Shape patch_size)
{
	int radius_x = patch_size.size_x / 2;
	int radius_y = patch_size.size_y / 2;

	Mask extended_domain = inpainting_domain;	// NOTE: implicit deep copy due to the ImmutableMask -> Mask casting

	FixedMask::iterator it;
	for (it = inpainting_domain.begin(); it != inpainting_domain.end(); ++it) {
		int x = it->x;
		int y = it->y;

		bool is_surrounded = inpainting_domain.test(x + 1, y) &&
							 inpainting_domain.test(x - 1, y) &&
							 inpainting_domain.test(x, y + 1) &&
							 inpainting_domain.test(x, y - 1);

		if (!is_surrounded) {
			int x_a = max(x - radius_x, 0);
			int x_b = min(x + radius_x, (int)inpainting_domain.get_size_x() - 1);
			int y_a = max(y - radius_y, 0);
			int y_b = min(y + radius_y, (int)inpainting_domain.get_size_y() - 1);

			for (int i = x_a; i <= x_b; i++) {
				for (int j = y_a; j <= y_b; j++) {
					if (!inpainting_domain.get(i, j)) {
						extended_domain.mask(i, j);
					}
				}
			}
		}
	}

	return extended_domain;
}


/**
 * Computes the confidence mask in such a way, that pixels outside the given mask have confidence of 1.0
 * and inside the masked region confidence values gradually decay up to the given asymptotic value.
 *
 * @param domain Masked region
 * @param decay_time Controls the speed of confidence values decay
 * @param asymptotic_value Lower boundary for confidence values
 */
Image<float> ImageInpainting::calculate_confidence_mask(FixedMask domain, float decay_time, float asymptotic_value)
{
	Image<float> confidence_mask(domain.get_size(), 1.0f);
	if (decay_time > 0) {
		Image<float> distances_to_boundary = DistanceTransform::calculate(domain);
		FixedMask::iterator it;
		for (it = domain.begin(); it != domain.end(); ++it) {
			float distance = distances_to_boundary(it->x, it->y);
			float value = ( 1 - asymptotic_value ) * exp( -distance / decay_time ) + asymptotic_value;
			confidence_mask(it->x, it->y) = value;
		}
	} else {
		FixedMask::iterator it;
		for (it = domain.begin(); it != domain.end(); ++it) {
			confidence_mask(it->x, it->y) = asymptotic_value;
		}
	}

	return confidence_mask;
}


/**
 * Adds a margin (unmasked points) of the given width at the border of two given masks.
 * Normally one mask should be a source region mask and another - target region mask (order does not matter).
 *
 * @param first_mask First mask to be processed
 * @param second_mask Second mask to be processed
 * @param margin Margin width to be added
 */
void ImageInpainting::add_margin(Mask first_mask, Mask second_mask, int margin)
{
	for (int i = 0; i < margin; i++) {
		for (uint x = 0; x < first_mask.get_size_x(); x++) {
			first_mask.unmask(x, i);
			first_mask.unmask(x, first_mask.get_size_y() - i - 1);
			second_mask.unmask(x, i);
			second_mask.unmask(x, second_mask.get_size_y() - i - 1);
		}
		for (uint y = 0; y < first_mask.get_size_y(); y++) {
			first_mask.unmask(i, y);
			first_mask.unmask(first_mask.get_size_x() - i - 1, y);
			second_mask.unmask(i, y);
			second_mask.unmask(second_mask.get_size_x() - i - 1, y);
		}
	}
}


inline int ImageInpainting::round(float value)
{
	return (value > 0.0) ? floor(value + 0.5) : ceil(value - 0.5);
}




