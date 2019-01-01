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

#include <string>
#include <ctime>

#include "image_inpainting.h"
#include "patch_match.h"
#include "a_patch_distance.h"
#include "patch_non_local_means.h"
#include "patch_non_local_medians.h"
#include "patch_non_local_poisson.h"
#include "patch_non_local_depth.h"
#include "l1_norm_patch_distance.h"
#include "l2_norm_patch_distance.h"
#include "l2_combined_patch_distance.h"
#include "io_utility.h"

#include <stdlib.h>
#include <stdio.h>
#include <stdexcept>

// NOTE: define DBG_OUTPUT in image_inpainting.h header to turn on the output after each iteration and other debug images output

using namespace std;

/**
 * @param c pointer to original argc
 * @param v pointer to original argv
 * @param name option name after hyphen
 * @param default_value default value (if NULL, the option takes no argument)
 */
static const char *pick_option(int *c, char ***v, const char *name, const char *default_value)
{
	int argc = *c;
	char **argv = *v;
	int id = default_value ? 1 : 0;
	for (int i = 0; i < argc - id; i++) {
		if (argv[i][0] == '-' && 0 == strcmp(argv[i] + 1, name)) {
			char *r = argv[i + id] + 1 - id;
			*c -= (id + 1);
			for (int j = i; j < argc - id; j++) {
				(*v)[j] = (*v)[j + id + 1];
			}

			return r;
		}
	}

	return default_value;
}



int main(int argc, char *argv[])
{
	// get parameters from command line
	int patch_side						= atoi(pick_option(&argc, &argv, "patch"  , "9"));			// 7
	int inpainting_iterations			= atoi(pick_option(&argc, &argv, "iters"  , "300"));		// 50 ok too
	string method_name					=      pick_option(&argc, &argv, "method" , "nlmeans");		// nlpoisson, nlmedians or nlmeans
	int scales_amount					= atoi(pick_option(&argc, &argv, "scales" , "7"));
	float coarsest_rate					= atof(pick_option(&argc, &argv, "coarse" , "0"));
	float confidence_decay_time			= atof(pick_option(&argc, &argv, "conft"  , "5.0"));
	float confidence_asymptotic_value	= atof(pick_option(&argc, &argv, "confa"  , "0.1"));
	float lambda						= atof(pick_option(&argc, &argv, "lambda" , "0.05"));
	string initialization_type_name		=      pick_option(&argc, &argv, "init"   , "poisson");		// poisson, avg, black, none
	float patch_sigma					= atof(pick_option(&argc, &argv, "psigma" , "10000.0"));	// uniform weights
	string show_nnf_file				=      pick_option(&argc, &argv, "shownnf", "");
	string show_pyramid_file			=      pick_option(&argc, &argv, "showpyr", "");
	string depth_name			=      pick_option(&argc, &argv, "depth", "");                  // Depth filename

	if (argc < 4) {
		// display usage message and quit
		fprintf(stderr, "USAGE:\n");
		fprintf(stderr, "%s input mask output [OPTIONS]\n\n", argv[0]);
		fprintf(stderr, "Available options are:\n");
		fprintf(stderr, " -method \tmethod name [nlmeans/nlmedians/nlpoisson] (%s)\n", method_name.c_str());
		fprintf(stderr, " -patch  \tpatch side (%d)\n", patch_side);
		fprintf(stderr, " -iters  \tinpainting iterations (%d)\n", inpainting_iterations);
		fprintf(stderr, " -scales \tscales amount (%d)\n", scales_amount);
		fprintf(stderr, " -coarse \tcoarsest rate (%g)\n", coarsest_rate);
		fprintf(stderr, " -conft  \tconfidence decay time (%g)\n", confidence_decay_time);
		fprintf(stderr, " -confa  \tconfidence asymptotic value (%g)\n", confidence_asymptotic_value);
		fprintf(stderr, " -lambda \tlambda (%g)\n", lambda);
		fprintf(stderr, " -init   \tinitialization type [poisson/black/avg/none] (%s)\n", initialization_type_name.c_str());
		fprintf(stderr, " -psigma \tGaussian patch weights (%g)\n", patch_sigma);
		fprintf(stderr, " -showpyr\tPREFIX write intermediate pyramid results\n");
		fprintf(stderr, " -shownnf\tFILENAME write illustration of the final NNF\n");
		return 1;
	}

	// get image paths from command line (or use default values)
	int i=1;
	string input_name  = (argc>i) ? argv[i]: "a057.png"; i++;
	string mask_name   = (argc>i) ? argv[i]: "a057_msk.png"; i++;
	string output_name = (argc>i) ? argv[i]: "a057_result"; i++;
	//printf("%s %s %s\n", input_name.c_str(), mask_name.c_str(), output_name.c_str());

#ifdef DBG_OUTPUT
	// define prefix for debug output
	IOUtility::set_prefix("dbg/");
#endif

	// read image and mask
	Image<float> input = IOUtility::rgb_to_lab(IOUtility::read_rgb_image(input_name));
	Image<float> mask_image = IOUtility::read_mono_image(mask_name);

	if (input.get_size() != mask_image.get_size()) {
		throw std::runtime_error("ERROR: Input image and mask should be of the same size");
	}

	// convert mask image to a Mask object
	bool is_mask_empty = true;
	Mask mask(mask_image.get_size_x(), mask_image.get_size_y());
	for (uint y = 0; y < mask_image.get_size_y(); y++) {
		for (uint x = 0; x < mask_image.get_size_x(); x++) {
			if (mask_image(x, y) > 0) {
				mask.mask(x, y);
				is_mask_empty = false;
			}
		}
	}

	if (is_mask_empty) {
		throw std::runtime_error("ERROR: Mask should not be empty (all black)");
	}

    // read depth map (if applicable)
	bool depth_exists = false;
	Image<float> depth;
	if (!depth_name.empty())    {
		depth_exists = true;
		depth = IOUtility::read_mono_image(depth_name);
		input = IOUtility::cat(input, depth);
	}
	// DEBUG: Depth concatenate test
	
	// if (depth_exists) {
	// 	Image<float> all;
	// 	Image<float> separated_RGB(input.get_size(), (uint)3);
	// 	all = IOUtility::cat(input, depth);
	// 	separated_RGB = IOUtility::lab_to_rgb(all);
	// 	IOUtility::write_rgb_image("test.png", separated_RGB);
	// }

	// rest of parameters

	// default value for the coarsest scale, using Alasdair Newson rule of thumb
	//
	//	  "... the size of the coarsest image has to be such that that the
	//	   patch size is no smaller than twice the size of the coarsest
	//	   inpainting domain ..." 
	//
	if (coarsest_rate == 0)
	{
		Image<float> distances_to_boundary = DistanceTransform::calculate(mask);

		float max_dist = 0;
		for (uint y = 0; y < distances_to_boundary.get_size_y(); y++)
			for (uint x = 0; x < distances_to_boundary.get_size_x(); x++)
				max_dist = (distances_to_boundary(x, y) > max_dist)
				         ?  distances_to_boundary(x, y) : max_dist;

		coarsest_rate = std::min(1.5f*patch_side/max_dist, 1.f);
	}

	// check errors in user given parameters
	if (patch_side <= 0) {
		throw std::runtime_error("ERROR: patch_side needs to be a positive integer.");
	}
	if ((coarsest_rate > 1) || (coarsest_rate < 0)) {
		throw std::runtime_error("ERROR: coarsest rate should be between 0 and 1.");
	}
	if (inpainting_iterations < 0) {
		throw std::runtime_error("ERROR: inpainting_iterations cannot be negative.");
	}
	if (confidence_decay_time <= 0.f) {
		throw std::runtime_error("ERROR: confidence_decay_time needs to be positive.");
	}
	if ((lambda > 1.f) || (lambda < 0.f)) {
		throw std::runtime_error("ERROR: lambda should be between 0 and 1.");
	}
	if (scales_amount < 1) {
		throw std::runtime_error("ERROR: number of scales must be at least 1.");
	}

	// if coarsest scale size ratio is equal to 1, set the number of scales to 1
	if (std::abs(coarsest_rate - 1.0f) < 0.00001f) {
		scales_amount = 1;
		printf("\tnote: since the coarsest scale size ratio is 1.0, the number of scales is set to 1\n");
	}

	// set initialization type
	ImageInpainting::InitType init_type;
	if (initialization_type_name.compare("poisson") == 0) {
		init_type = ImageInpainting::InitPoisson;
    } else if (initialization_type_name.compare("black") == 0) {
		init_type = ImageInpainting::InitBlack;
	} else if (initialization_type_name.compare("avg") == 0) {
		init_type = ImageInpainting::InitAvg;
	} else if (initialization_type_name.compare("none") == 0) {
		init_type = ImageInpainting::InitNone;
	} else {
		throw std::runtime_error("ERROR: Unknown initialization type");
	}

	// define inpainting parameters
	float tolerance = 0.1;
	float subsampling_rate = ImageInpainting::calculate_subsampling_rate(coarsest_rate, scales_amount);

	// define patch distance calculation parameters
	Shape weights_update_patch_size = Shape(patch_side, patch_side);
	float weights_update_sigma = patch_sigma;	

	// define Patch-Match parameters
	int patch_match_iterations = 10;
	int random_shots_limit = 30;

	// define image updating parameters
	Shape image_update_patch_size = Shape(patch_side, patch_side);
	float image_update_sigma = weights_update_sigma;

	// initialize multiscale inpainting algorithm
	ImageInpainting image_inpainting = ImageInpainting(inpainting_iterations,
	                                                   tolerance,
	                                                   scales_amount,
	                                                   subsampling_rate,
	                                                   confidence_decay_time,
	                                                   confidence_asymptotic_value,
	                                                   init_type);

	// create PatchDistance and ImageUpdating objects
	APatchDistance *patch_distance;
	AImageUpdating *image_updating;
	if(method_name.compare("nlmeans") == 0) {
		image_updating = new PatchNonLocalMeans(image_update_patch_size, image_update_sigma);
		patch_distance = new L2NormPatchDistance(weights_update_patch_size, weights_update_sigma);
	} else if(method_name.compare("nlmedians") == 0) {
		image_updating = new PatchNonLocalMedians(image_update_patch_size, image_update_sigma);
		patch_distance = new L1NormPatchDistance(weights_update_patch_size, weights_update_sigma);
	} else if (method_name.compare("nlpoisson") == 0) {
		image_updating = new PatchNonLocalPoisson(image_update_patch_size, image_update_sigma, lambda, 0.000001, 1000);	// 0.000001, 1000
		patch_distance = new L2CombinedPatchDistance(lambda, weights_update_patch_size, weights_update_sigma);
	} else if (method_name.compare("nldepth") == 0) {
		// TODO: image update method for RGBD
		image_updating = new PatchNonLocalDepth(image_update_patch_size, image_update_sigma, lambda, 0.000001, 1000);	// 0.000001, 1000
		patch_distance = new L2CombinedPatchDistance(lambda, weights_update_patch_size, weights_update_sigma);
	} else {
		throw std::runtime_error("ERROR: Unknown method name");
	}

	// create PatchMatch object
	PatchMatch *patch_match = new PatchMatch(patch_distance, patch_match_iterations, random_shots_limit, -1);

	// link PacthMatch and ImageUpdating objects to multiscale image inpainter
	image_inpainting.set_weights_updating(patch_match);
	image_inpainting.set_image_updating(image_updating);

	// init random generator (for PatchMatch)
	srand(time(NULL));
	
	// tell algorithm to keep original image pyramid and nnf pyramid, if needed
	image_inpainting.keep_intermediate(!show_nnf_file.empty() || !show_pyramid_file.empty());

	// do multiscale inpainting
	Image<float> output = image_inpainting.process(input, mask);

	// clean
	delete patch_match;
	delete patch_distance;
	delete image_updating;

	// compose output name
	stringstream output_str;
	output_str << output_name;
	output_str << "_" << method_name;
	if (method_name == "nlpoisson") {
		output_str << "_l" << lambda;
	}
	output_str << "_sc" << scales_amount << "_" << coarsest_rate;
	if (init_type == ImageInpainting::InitAvg) {
		output_str << "_initavg";
	} else if (init_type == ImageInpainting::InitBlack) {
		output_str << "_init000";
	} else {
		output_str << "_initnone";
	}
	output_str << "_ps" << patch_side << "_" << weights_update_sigma;
	output_str << "_conf" << confidence_decay_time << "_" << confidence_asymptotic_value;
	output_str << ".png";

	// save result
#ifndef IPOL_DEMO
	IOUtility::write_rgb_image(output_str.str(), IOUtility::lab_to_rgb(output));
	if (depth_exists)
		IOUtility::write_rgb_image("inpainted_depth.png", IOUtility::get_depth(output));
#else
   // in the IPOL demo the output filename is fixed just output_name
	IOUtility::write_rgb_image(output_name, IOUtility::lab_to_rgb(output));
#endif
   
	// save the output pyramid
	if (!show_pyramid_file.empty()) {
		vector<Image<float > > output_pyramid = image_inpainting.get_output_image_pyramid();
		for (uint n=0; n < output_pyramid.size(); n++) {
			stringstream output_s;
			output_s << show_pyramid_file;
			output_s << n;
			output_s << ".png";			
			IOUtility::write_rgb_image(output_s.str(), IOUtility::lab_to_rgb(output_pyramid[n]));
		}
	}

	
	// save a representation of the NNF
	if (!show_nnf_file.empty()) {
		vector<Image<Point> > nnf_pyramid = image_inpainting.get_nnf_pyramid();
      Image<Point> lastNNF = nnf_pyramid[nnf_pyramid.size()-1];
		Image<float> show_nnf = IOUtility::lab_to_rgb(output);

		for (uint c = 0; c < show_nnf.get_number_of_channels(); c++) {
			for (uint y = 1; y < show_nnf.get_size_y(); y++) {
				for (uint x = 1; x < show_nnf.get_size_x(); x++) {

               // show as boundary 
					if ( mask(x,y) != mask(x-1,y) || mask(x,y) != mask(x,y-1) ) {
						show_nnf(x,y,c) = 0;
					}
					if (mask(x,y)) {
						if (x > 0) {
                     Point p2(x,y), p1(x-1,y);
							if (lastNNF(x-1,y) - p1 != lastNNF(x,y) - p2 ) {
								show_nnf(x,y,c) = 0;
							}
						}
						if (y > 0) {
                     Point p2(x,y), p1(x,y-1);
							if (lastNNF(x,y-1) - p1 != lastNNF(x,y) - p2 ) {
								show_nnf(x,y,c) = 0;
							}
						}
					}

               // show as color displacement map
//					if (mask(x,y)) {
//						float lx = (float(lastNNF(x,y).x)/show_nnf.get_size_x());
//						float ly = (float(lastNNF(x,y).y)/show_nnf.get_size_y());
//
//						if(c==0)
//							show_nnf(x,y,c) = lx*255.0;
//						if(c==1)
//							show_nnf(x,y,c) = lx*ly*255.0;
//						if(c==2)
//							show_nnf(x,y,c) = ly*255.0;
//					}

				}
			}
		}
		IOUtility::write_rgb_image(show_nnf_file, show_nnf);
	}
}
