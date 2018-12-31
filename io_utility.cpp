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

#include "io_utility.h"

extern "C" {
#include "iio.h"
}


string IOUtility::_prefix = "";


/**
 * Reads image in PGM format from file with provided path\name
 * to the object of 'Image' class
 */
Image<float> IOUtility::read_pgm_image(const string &name)
{
	// open file
	FILE *f = fopen(name.data(),"rb");	// COMPATIBILITY: for win 'rb' file mode instead of just 'r'
	if( f == NULL ) {
		return Image<float>();
	}

	// read header
	bool isBinary = false;
	int c, x_size,y_size,depth;

	if ( fgetc(f) != 'P' ) {
		return Image<float>();
	}

	if( (c=fgetc(f)) == '2' ) {
		isBinary = false;
	} else if ( c == '5' ) {
		isBinary = true;
	} else {
		return Image<float>();
	}

	skip_spaces_and_comments(f);
	fscanf(f,"%d",&x_size);
	skip_spaces_and_comments(f);
	fscanf(f,"%d",&y_size);
	skip_spaces_and_comments(f);
	fscanf(f,"%d",&depth);

	// allocate memory
	Image<float> image(x_size, y_size);

	// read data
	skip_spaces_and_comments(f);
	for(int y = 0; y < y_size; y++) {
		for(int x = 0; x < x_size; x++) {
			int value = isBinary ? fgetc(f) : get_number(f);
			image(x, y) = value;
		}
	}

	// close and return
	fclose(f);

	return image;
}

/**
 * Writes image in PGM format to file with provided path\name
 * from the object of 'Image' class
 */
void IOUtility::write_pgm_image(const string &name, FixedImage<float> image)
{
	// open file
	FILE *f = fopen(name.data(),"wb");

	// write header
	putc('P', f);
	putc('5', f);

	// write attributes
	int size_x = image.get_size_x();
	int size_y = image.get_size_y();
	fprintf(f, "\n%d %d\n%d\n", size_x, size_y, 255);

	// write data
	if (image.get_number_of_channels() == 1) {
		for (int y = 0; y < size_y; y++) {
			for (int x = 0; x < size_x; x++) {
				char value = (char)image(x,y);
				putc(value, f);
			}
		}
	} else if (image.get_number_of_channels() == 3) {
		for (int y = 0; y < size_y; y++) {
			for (int x = 0; x < size_x;x++) {
				float red = image(x, y, 0);
				float green = image(x, y, 1);
				float blue = image(x, y, 2);
				// luminance value
				char value = (char)(0.2126 * red + 0.7152 * green + 0.0722 * blue);
				putc(value, f);
			}
		}
	}

	// close file
	fclose(f);
}


Image<float> IOUtility::read_mono_image(const string &name)
{
	int width, height;

	float *image_data = iio_read_image_float(name.data(), &width, &height);

	if (! image_data) {
		fprintf(stderr, "read_mono_image: cannot read %s\n", name.data());
		exit(1);
	}

	Image<float> image(width, height);

	// copy from image_data to Image<float>
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int index = (y * width + x);
			float value = image_data[index];
			image(x, y) = value;
		}
	}

	free(image_data);

	return image;
}


Image<float> IOUtility::read_rgb_image(const string &name)
{
	int width, height, channels;

	float *image_data = iio_read_image_float_vec(name.data(), &width, &height, &channels);

   if (! image_data) {
      fprintf(stderr, "read_rgb_image: cannot read %s\n", name.data());
      exit(1);
   } 

	Image<float> image(width, height, (uint)3);

	// copy from image_data to Image<float>
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int index = channels * (y * width + x);
			image(x, y, 0) = image_data[index];
			image(x, y, 1) = channels>1 ? image_data[index + 1]: image_data[index];
			image(x, y, 2) = channels>2 ? image_data[index + 2]: image_data[index];
		}
	}

	free(image_data);

	return image;
}


void IOUtility::write_rgb_image(const string &name, FixedImage<float> image)
{
	int width = image.get_size_x();
	int height = image.get_size_y();

	float *image_data = new float[width * height * 3];

	// copy from Image<float> to image_data
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int index = 3 * (y * width + x);
			image_data[index + 0] = image(x, y, 0);
			image_data[index + 1] = image(x, y, 1);
			image_data[index + 2] = image(x, y, 2);
		}
	}

	iio_save_image_float_vec(const_cast<char* >(name.data()), image_data, width, height, 3);

   delete[] image_data;
}

void IOUtility::write_mono_image(const string &name, FixedImage<float> image)
{
	int width = image.get_size_x();
	int height = image.get_size_y();

	float *image_data = new float[width * height];

	// copy from ImmutableImage<float> to image_data
	if (image.get_number_of_channels() == 1) {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int index = (y * width + x);
				image_data[index] = image(x, y);
			}
		}
	} else if (image.get_number_of_channels() == 3) {
		for (int y = 0; y < height; y++) {
			for (int x = 0; x < width; x++) {
				int index = (y * width + x);
				float red = image(x, y, 0);
				float green = image(x, y, 1);
				float blue = image(x, y, 2);
				image_data[index] = 0.2126 * red + 0.7152 * green + 0.0722 * blue;
			}
		}
	}

	iio_save_image_float(const_cast<char* >(name.data()), image_data, width, height);
	delete[] image_data;
}


/*
 * Sets all masked points to white (255) and unmasked - to black (0)
 */
Image<float> IOUtility::mask_to_greyscale(FixedMask mask)
{
	Shape size = mask.get_size();
	Image<float> result(size.size_x, size.size_y, 0.0f);

	FixedMask::iterator it;
	for (it = mask.begin(); it != mask.end(); ++it) {
		result(it->x, it->y) = 255.0;
	}

	return result;
}


/*
 * Translates range [0.0; 1.0] into the range [0; 255]
 * Might be used to draw probabilities or 'soft' masks
 */
Image<float> IOUtility::probability_to_greyscale(FixedImage<float> probabilities)
{
	Shape size = probabilities.get_size();
	Image<float> result(size.size_x, size.size_y, 0.0f);

	for (unsigned int y = 0; y < size.size_y; y++) {
		for (unsigned int x = 0; x < size.size_x; x++) {
			float value = fmax(0.0, fmin(1.0, probabilities(x, y)));
			result(x, y) = 255.0 * value;
		}
	}

	return result;
}


Image<float> IOUtility::rgb_to_lab(FixedImage<float> image)
{
	if (image.get_number_of_channels() != 3) {
		return image;
	}

	Image<float> image_lab(image.get_size(), 3, 0.0f);
	const float* data_rgb = image.raw();
	float* data_lab = image_lab.raw();
	long number_of_pixels = image.get_size_x() * image.get_size_y();
	float *xyz_buffer = new float[3];

	//#pragma omp parallel for
	for (int i = 0; i < number_of_pixels; i++) {
		rgb_to_xyz(data_rgb + i * 3, xyz_buffer);
		xyz_to_lab(xyz_buffer, data_lab + i * 3);
	}

	delete[] xyz_buffer;

	return image_lab;
}


Image<float> IOUtility::lab_to_rgb(FixedImage<float> image)
{
	if (image.get_number_of_channels() != 3 && image.get_number_of_channels() != 4) {
		return image;
	}
    
	Image<float> image_rgb(image.get_size(), 3, 0.0f);
	const float* data_lab = image.raw();
	float* data_rgb = image_rgb.raw();
	long number_of_pixels = image.get_size_x() * image.get_size_y();
	float *xyz_buffer = new float[3];
    if (image.get_number_of_channels() == 4)    {
		Image<float> img(image.get_size(), (uint)3);
		Image<float> depth(image.get_size(), (uint)1);
		IOUtility::separate(image, img, depth);
		data_lab = img.raw();
	} 
	//#pragma omp parallel for
	for (int i = 0; i < number_of_pixels; i++) {
		lab_to_xyz(data_lab + i * 3, xyz_buffer);
		xyz_to_rgb(xyz_buffer, data_rgb + i * 3);
	}

	delete[] xyz_buffer;

	return image_rgb;
}


string IOUtility::compose_file_name(const string &name)
{
	return _prefix + name;
}


string IOUtility::compose_file_name(const string &name, int index, const string &extension)
{
	stringstream stream;
	stream << _prefix << name << "_" << index << "." << extension;
	string file_name = stream.str();
	return file_name;
}


string IOUtility::compose_file_name(const string &name, int index, int index2, const string &extension)
{
	stringstream stream;
	stream << _prefix << name << "_" << index << "_" << setw(3) << setfill('0') << index2 << "." << extension;
	string file_name = stream.str();
	return file_name;
}


void IOUtility::set_prefix(const string &prefix)
{
	_prefix = prefix;
}

/* Private */

void IOUtility::skip_spaces_and_comments(FILE * f)
{
  int c;
  do
    {
      while(isspace(c=fgetc(f))); // skip spaces
      if(c=='#') while((c=fgetc(f))!='\n'); // skip comments
    }
  while(c == '#');
  ungetc(c,f);
}

/**
 *  Reads a number digit by digit.
 */
int IOUtility::get_number(FILE * f)
{
  int num, c;

  while(isspace(c=fgetc(f)));
  if(!isdigit(c)) exit(1);
  num = c - '0';
  while( isdigit(c=fgetc(f)) ) num = 10 * num + c - '0';

  return num;
}


void IOUtility::rgb_to_xyz(const float *rgb, float *xyz)
{
	float aux_r = rgb[0] / 255.0f;
	float aux_g = rgb[1] / 255.0f;
	float aux_b = rgb[2] / 255.0f;

	aux_r = (aux_r > 0.04045) ? pow((aux_r + 0.055f) / 1.055f , 2.4f) : aux_r / 12.92;
	aux_g = (aux_g > 0.04045) ? pow((aux_g + 0.055f) / 1.055f , 2.4f) : aux_g / 12.92;
	aux_b = (aux_b > 0.04045) ? pow((aux_b + 0.055f) / 1.055f , 2.4f) : aux_b / 12.92;

	aux_r *= 100.0f;
	aux_g *= 100.0f;
	aux_b *= 100.0f;

	xyz[0] = aux_r * 0.412453f + aux_g * 0.357580f + aux_b * 0.180423f;
	xyz[1] = aux_r * 0.212671f + aux_g * 0.715160f + aux_b * 0.072169f;
	xyz[2] = aux_r * 0.019334f + aux_g * 0.119193f + aux_b * 0.950227f;
}


void IOUtility::xyz_to_lab(const float *xyz, float *lab)
{
	// normalize by the reference white
	float aux_x = xyz[0] / 95.047f;
	float aux_y = xyz[1] / 100.000f;
	float aux_z = xyz[2] / 108.883f;

	aux_x = (aux_x > 0.008856f) ? pow(aux_x, 1.0f / 3.0f) : (7.787f * aux_x) + (16.0f / 116.0f);
	aux_y = (aux_y > 0.008856f) ? pow(aux_y, 1.0f / 3.0f) : (7.787f * aux_y) + (16.0f / 116.0f);
	aux_z = (aux_z > 0.008856f) ? pow(aux_z, 1.0f / 3.0f) : (7.787f * aux_z) + (16.0f / 116.0f);

	lab[0] = (116.0f * aux_y) - 16.0f;
	lab[1] = 500.0f * (aux_x - aux_y);
	lab[2] = 200.0f * (aux_y - aux_z);
}


void IOUtility::lab_to_xyz(const float *lab, float *xyz)
{
	float aux_y = (lab[0] + 16.0f) / 116.0f;
	float aux_x = lab[1] / 500.0 + aux_y;
	float aux_z = aux_y - lab[2] / 200.0f;

	aux_x = (pow(aux_x, 3.0f) > 0.008856f) ? pow(aux_x, 3.0f) : (aux_x - 16.0f / 116.0f) / 7.787f;
	aux_y = (pow(aux_y, 3.0f) > 0.008856f) ? pow(aux_y, 3.0f) : (aux_y - 16.0f / 116.0f) / 7.787f;
	aux_z = (pow(aux_z, 3.0f) > 0.008856f) ? pow(aux_z, 3.0f) : (aux_z - 16.0f / 116.0f) / 7.787f;

	xyz[0] = aux_x * 95.047f;
	xyz[1] = aux_y * 100.000f;
	xyz[2] = aux_z * 108.883f;
}


void IOUtility::xyz_to_rgb(const float *xyz, float *rgb)
{
	float aux_x = xyz[0] / 100.0f;
	float aux_y = xyz[1] / 100.0f;
	float aux_z = xyz[2] / 100.0f;

	float aux_r = aux_x *  3.240479f + aux_y * -1.537150f + aux_z * -0.498535f;
	float aux_g = aux_x * -0.969256f + aux_y *  1.875992f + aux_z *  0.041556f;
	float aux_b = aux_x *  0.055648f + aux_y * -0.204043f + aux_z *  1.057311f;

	aux_r = (aux_r > 0.0031308f) ? 1.055f * pow(aux_r, 1.0f / 2.4f) - 0.055f : 12.92f * aux_r;
	aux_g = (aux_g > 0.0031308f) ? 1.055f * pow(aux_g, 1.0f / 2.4f) - 0.055f : 12.92f * aux_g;
	aux_b = (aux_b > 0.0031308f) ? 1.055f * pow(aux_b, 1.0f / 2.4f) - 0.055f : 12.92f * aux_b;

	rgb[0] = aux_r * 255.0f;
	rgb[1] = aux_g * 255.0f;
	rgb[2] = aux_b * 255.0f;
}

// Concatenate RGB and Depth Channel
Image<float> IOUtility::cat(FixedImage<float> image, FixedImage<float> depth)    {
	Image<float> out(image.get_size(), (uint)4);
	if(image.get_number_of_channels() != 3 && depth.get_number_of_channels() != 1)    {
		throw std::runtime_error("cat: Channels do not match.");
	}
	if(image.get_size() != depth.get_size())    {
		throw std::runtime_error("cat: Size does not match.");
	}
	for (uint row = 0; row < image.get_size_y(); row++)    {
		for (uint col = 0; col < image.get_size_x(); col++)    {
			uint ch = 0;
			for (ch = 0; ch < 3; ch++)    {
				out.at(col, row, ch) = image.at(col, row, ch);
			}
			out.at(col, row, ch) = depth.at(col, row);
		}
	}
	return out;
}
// Separate depth and RGB
void IOUtility::separate(FixedImage<float> input, Image<float> &rgb, Image<float> &depth)    {
	if (input.get_number_of_channels() != 4)    {
		throw std::runtime_error("Separate: Channels do not match.");
	}
	for (uint row = 0; row < input.get_size_y(); row++)    {
		for (uint col = 0; col < input.get_size_x(); col++)    {
			uint ch = 0;
			for (ch = 0; ch < 3; ch++)    {
				rgb.at(col, row, ch) = input.at(col, row, ch);
			}
			depth.at(col, row) = input.at(col, row, ch);
		}
	}
}

Image<float> IOUtility::get_depth(FixedImage<float> input)    {
	if (input.get_number_of_channels() != 4)    {
		throw std::runtime_error("Separate: Channels do not match.");
	}
	Image<float> rgb(input.get_size(), (uint)3);
	Image<float> depth(input.get_size());
	IOUtility::separate(input, rgb, depth);
	return depth;
}