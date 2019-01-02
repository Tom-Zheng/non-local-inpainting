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

#ifndef IO_UTILITY_H_
#define IO_UTILITY_H_

#include <string>
#include <sstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

#include "image.h"
#include "mask.h"

/**
 * Contains methods for reading and writings images as well as
 * some auxiliary methods. Acts as a proxy for IIO library.
 */
class IOUtility
{
public:
	// Reading and writing of a grayscale image without IIO.
	static Image<float> read_pgm_image(const string &name);
	static void write_pgm_image(const string &name, FixedImage<float> image);

	// Reading and writing of a grayscale image using IIO.
	static Image<float> read_mono_image(const string &name);
	static void write_mono_image(const string &name, FixedImage<float> image);

	// Reading and writing for a colored image using IIO.
	static Image<float> read_rgb_image(const string &name);
	static void write_rgb_image(const string &name, FixedImage<float> image);

	// Sets all masked points to white (255) and unmasked - to black (0)
	static Image<float> mask_to_greyscale(FixedMask mask);

	// Translates range [0.0; 1.0] into the range [0; 255]
	static Image<float> probability_to_greyscale(FixedImage<float> probabilities);

	static Image<float> rgb_to_lab(FixedImage<float> image);
	static Image<float> lab_to_rgb(FixedImage<float> image);

	static string compose_file_name(const string &name);
	static string compose_file_name(const string &name, int index, const string &extension);
	static string compose_file_name(const string &name, int index, int index2, const string &extension);

	// Prefix is applied to a file name by compose_file_name functions
	static void set_prefix(const string &prefix);

	// Concatenate RGB and Depth Channel
	static Image<float> cat(FixedImage<float> image, FixedImage<float> depth);
	static void cat(FixedImage<float> image, FixedImage<float> depth, Image<float> out);
	// Separate depth and RGB
	static void separate(FixedImage<float> input, Image<float> &rgb, Image<float> &depth);
	static Image<float> get_depth(FixedImage<float> input);

private:
	static string _prefix;

	static void skip_spaces_and_comments(FILE * f);
	static int get_number(FILE * f);

	static void rgb_to_xyz(const float *rgb, float *xyz);
	static void xyz_to_lab(const float *xyz, float *lab);
	static void lab_to_xyz(const float *lab, float *xyz);
	static void xyz_to_rgb(const float *xyz, float *rgb);

};




#endif /* IO_UTILITY_H_ */
