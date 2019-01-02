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

#ifndef IMAGE_INPAINTING_H_
#define IMAGE_INPAINTING_H_

#include <vector>
#include <limits>
#include <stdio.h>

#include "image.h"
#include "mask.h"
#include "patch_match.h"
#include "shape.h"
#include "sampling.h"
#include "distance_transform.h"
#include "a_image_updating.h"

/// define DBG_OUTPUT to turn on the output after each iteration and other debug images output
// #define DBG_OUTPUT

// #ifdef DBG_OUTPUT
// #include "io_utility.h"
// #endif

/**
 * The main class implementing image inpainting method.
 */
class ImageInpainting
{
public:
	enum InitType { InitBlack, InitAvg, InitNone, InitPoisson };

	ImageInpainting();
	ImageInpainting(int iterations_amount,
					float tolerance,
					int scales_amount,
					float subsampling_rate,
					float confidence_decay_time,
					float confidence_asymptotic_value,
					InitType initialization_type);

	static float calculate_subsampling_rate(float size_ratio, unsigned int scales_amount);

	Image<float> process(FixedImage<float> in, FixedMask mask);

	/// getters and setters for parameters
	int get_iterations_amount();
	void set_iterations_amount(int amount);
	float get_tolerance();
	void set_tolerance(float value);
	int get_scales_amount();
	void set_scales_amount(int amount);
	float get_subsampling_rate();
	void set_subsampling_rate(float rate);
	float get_confidence_decay_time();
	void set_confidence_decay_time(float value);
	float get_confidence_asymptotic_value();
	void set_confidence_asymptotic_value(float value);
	PatchMatch* get_weights_updating();
	void set_weights_updating(PatchMatch *patch_match);
	AImageUpdating* get_image_updating();
	void set_image_updating(AImageUpdating *image_updating);

	void keep_intermediate(bool value = true);
	vector<Image<float> > get_input_image_pyramid();
	vector<Image<float> > get_output_image_pyramid();
	vector<Mask> get_mask_pyramid();
	vector<Image<Point> > get_nnf_pyramid();



#ifdef DBG_OUTPUT
	int cm_ind;
#endif

private:
	static const float MASK_SAMPLING_THRESHOLD;

	// PatchMatch algorithm
	PatchMatch *_patch_match;

	// image updater
	AImageUpdating *_image_updating;

	// number of iterations
	int _iterations_amount;
	float _tolerance;

	// multiscale parameters
	int _scales_amount;      // number of scales
	float _subsampling_rate; // subsamplingn rate between scales

	// confidence weights parameters
	float _confidence_decay_time;
	float _confidence_asymptotic_value;

	// coarsest scale initialization (average, black or none)
	InitType _initialization_type;

	bool _keep_intermediate;
	vector<Image<float> > _original_image_pyramid;
	vector<Image<float> > _image_pyramid;
	vector<Mask> _mask_pyramid;
	vector<Image<Point> > _nnf_pyramid;

	// inpainting of one scale
	void inpaint_internal(Image<float> image,
						  FixedMask inpainting_domain,
						  FixedImage<float> confidence_mask,
						  FixedImage<Point> initial_nnf,
						  float tolerance);

	// sets all pixels in mask to color
	void initialize_with_color(Image<float> image,
							   FixedMask mask,
							   float* color);

	// returns average of known pixels OR black
	float* get_initial_color(Image<float> image,
							 FixedMask mask,
							 bool average);

	// extended domain are the centers of unknown patches
	Mask get_extended_domain(FixedMask inpainting_domain,
							 Shape patch_size);

	// upscaling from coarse scale to fine scale
	Image<Point> propagate_weights(Image<float> upper_level,
						   FixedMask upper_inpainting_domain,
						   FixedImage<float> upper_confidence_mask,
						   FixedImage<float> lower_level,
						   FixedMask lower_inpainting_domain);

	// computes confidence mask
	Image<float> calculate_confidence_mask(FixedMask domain,
										   float decay_time,
										   float asymptotic_value);

	// adds a margin to the image of half a patch
	void add_margin(Mask first_mask, Mask second_mask, int margin);

	// TODO: remove from here
	inline int round(float value);

};


#endif /* IMAGE_INPAINTING_H_ */
