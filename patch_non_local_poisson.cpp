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

#include "patch_non_local_poisson.h"

PatchNonLocalPoisson::PatchNonLocalPoisson()
	: AImageUpdating()
{
	_lambda = 0.5;
	_conjugate_gradient_tolerance = 0.0000001;
	_conjugate_gradient_iterations_limit = 1000;
}


PatchNonLocalPoisson::PatchNonLocalPoisson(Shape patch_size,
										   float gaussian_sigma,
										   float lambda,
										   float conjugate_gradient_tolerance,
										   int conjugate_gradient_iterations_limit)
	: AImageUpdating(patch_size, gaussian_sigma)
{
	_lambda = lambda;
	_conjugate_gradient_tolerance = conjugate_gradient_tolerance;
	_conjugate_gradient_iterations_limit = conjugate_gradient_iterations_limit;
}

double PatchNonLocalPoisson::update(Image<float> image,
									Image<float> original_image,
									FixedMask inpainting_domain,
									FixedMask extended_inpainting_domain,
									FixedImage<Point> nnf,
									FixedImage<float> confidence_mask)
{
	int size_x = image.get_size_x();
	int number_of_channels = image.get_number_of_channels();

	// calculate gradient.
	FixedImage<float> gradient = Gradient::calculate(image);

	// allocate memory
	int pixels_amount = image.get_size_x() * image.get_size_y();
	double *a1 = new double[pixels_amount];
	double *a2 = new double[pixels_amount];
	double **f1 = new double*[number_of_channels];
	double **f2 = new double*[number_of_channels];
	double **perturbation = new double*[number_of_channels];
	for (int i = 0; i < number_of_channels; i++) {
		f1[i] = new double[pixels_amount];
		f2[i] = new double[pixels_amount];
		perturbation[i] = new double[pixels_amount]();	// NOTE: default-initialized to 0.0
	}

	calculate_pde_coefficients(image, gradient, extended_inpainting_domain, nnf, confidence_mask, a1, a2, f1, f2);

	double** image_channels = split_image_into_channels(original_image);

	// get points from the inpainting domain
	vector<Point> points = inpainting_domain.get_masked_points();

	// solve Poisson PDE with Dirichlet boundary conditions using conjugate gradient
	for(int ch = 0; ch < number_of_channels; ch++) {	// iterate through the channels
		double* channel = new double[pixels_amount]();
		memcpy(channel, image_channels[ch], pixels_amount * sizeof(double));

		calculate_anisotropic_laplacian(channel, a2, inpainting_domain);

		vector<Point>::iterator it;
		for (it = points.begin(); it != points.end(); ++it) {
			int index = size_x * it->y + it->x;
			f1[ch][index] += - f2[ch][index] - a1[index] * image_channels[ch][index] + channel[index];

			perturbation[ch][index] = image(*it, ch) - image_channels[ch][index];
		}

		conjugte_gradient(perturbation[ch], a1, a2, f1[ch], inpainting_domain);

		delete[] channel;
	}

	// update the image
	double total_difference = 0.0;
	vector<Point>::iterator it;
	for (it = points.begin(); it != points.end(); ++it) {
		int index = size_x * it->y + it->x;

		for (int ch = 0; ch < number_of_channels; ch++) {
			float color_value = image_channels[ch][index] + perturbation[ch][index];

			// add to total difference
			float prev_value = image(*it, ch);
			total_difference += (prev_value - color_value) * (prev_value - color_value);

			image(*it, ch) = color_value;
		}
	}

	// free memory
	for (int i = 0; i < number_of_channels; i++) {
		delete[] f1[i];
		delete[] f2[i];
		delete[] perturbation[i];
		delete[] image_channels[i];
	}
	delete[] a1;
	delete[] a2;
	delete[] f1;
	delete[] f2;
	delete[] perturbation;
	delete[] image_channels;

	return total_difference;
}


/**
 * This function computes the coefficients of the following functional:
 *
 *    E(u,w) = \lambda*\int_{\widehat O} \int_{\Omega \setminus \widehat O}
 *             w(x,x') \int_{\Omega_p} g_a(y)(u(x + y) - u(x' + y) )^2 dydx'dx
 *             + (1-lambda)*\int_{\widehat O} \int_{\Omega \setminus \widehat O}
 *             w(x,x') \int_{\Omega_p} g_a(y)\|\nabla u(x + y) - \nabla u(x' + y) \|^2 dydx'dx
 * 	    + H(w)
 *
 * This function follows strictly the variational setting, thus:
 * the coefficients are computed based on a linear combination of the distance between
 * patches and the distance between patches of the gradient function.
 *
 * The patch distance is weighted by the intra-patch weight function.
 */
inline void PatchNonLocalPoisson::calculate_pde_coefficients(const FixedImage<float> &image,
														     const FixedImage<float> &gradient,
														     const FixedMask &inpainting_domain,
														     const FixedImage<Point> &nnf,
														     const FixedImage<float> &confidence_mask,
														     double *a1,
														     double *a2,
														     double **f1,
														     double **f2)
{
	int pixels_amount = image.get_size_x() * image.get_size_y();
	int number_of_channels = image.get_number_of_channels();
	Shape image_shape = image.get_size();

	// allocate auxiliary memory and set coefficients to zero
	memset(a1, 0, pixels_amount * sizeof(double));
	memset(a2, 0, pixels_amount * sizeof(double));
	double **F1 = new double*[number_of_channels];
	double **F2 = new double*[number_of_channels];
	for (int ch = 0; ch < number_of_channels; ch++) {
		F1[ch] = new double[pixels_amount]();
		F2[ch] = new double[pixels_amount]();

		memset(f1[ch], 0, pixels_amount * sizeof(double));
		memset(f2[ch], 0, pixels_amount * sizeof(double));
	}

	// calculate weights, if absent
	if (_patch_weighting.is_empty()) {
		_patch_weighting = GaussianWeights::calculate(_patch_size.size_x,
													  _patch_size.size_y,
													  _gaussian_sigma,
													  _gaussian_sigma);

	}

	int half_patch_size_x = _patch_size.size_x / 2;
	int half_patch_size_y = _patch_size.size_y / 2;

	// iterate through the inpainting domain
	FixedMask::iterator it;
	for (it = inpainting_domain.begin(); it != inpainting_domain.end(); ++it) {
		// restrict current patch to the image domain
		int x_first = max(it->x - half_patch_size_x, 0);
		int x_last  = min(it->x + half_patch_size_x, (int)image_shape.size_x - 1);
		int y_first = max(it->y - half_patch_size_y, 0);
		int y_last  = min(it->y + half_patch_size_y, (int)image_shape.size_y - 1);

		Point neighbor = nnf(*it);

		// NOTE: outside the inpainting domain Confidence is 1.0
		float confidence = (confidence_mask.is_not_empty() && inpainting_domain.test(*it)) ?
				confidence_mask(*it) : 1.0;

		for (int x = x_first; x <= x_last; x++) {
			for (int y = y_first; y <= y_last; y++) {
				// current point inside a patch centered at neighbor
				int n_x = neighbor.x + x - it->x;
				int n_y = neighbor.y + y - it->y;

				if (image_shape.contains(n_x, n_y)) {
					int index = image_shape.size_x * y + x;
					double weight = confidence * _patch_weighting(x - it->x + half_patch_size_x, y - it->y + half_patch_size_y);

					a1[index] += weight *      _lambda;
					a2[index] += weight * (1 - _lambda);

					for (int ch = 0; ch < number_of_channels; ch ++) {
						f1[ch][index] += weight * image(n_x, n_y, ch) * _lambda;

						F1[ch][index] += weight * gradient(n_x, n_y, ch * 2    ) * (1 - _lambda);
						F2[ch][index] += weight * gradient(n_x, n_y, ch * 2 + 1) * (1 - _lambda);
					}
				}
			}
		}

	}

	// compute divergence of the field F1,F2
	if (_lambda != 1) {
		for (int i = 0; i < number_of_channels; i++) {
			calculate_divergence(F1[i], F2[i], inpainting_domain, false, false, f2[i]);
		}
	}

	// free memory
	for (int i = 0; i < number_of_channels; i++) {
		delete[] F1[i];
		delete[] F2[i];
	}
	delete[] F1;
	delete[] F2;
}


inline void PatchNonLocalPoisson::calculate_divergence(const double *field_x,
													   const double *field_y,
													   const FixedMask &mask,
													   bool is_x_forward,
													   bool is_y_forward,
													   double *divergence)
{
	int size_x = mask.get_size_x();
	int size_y = mask.get_size_y();

	FixedMask::iterator it;
	for (it = mask.begin(); it != mask.end(); ++it) {
		int index = size_x * it->y + it->x;

		if (is_x_forward) {
			if (it->x == 0) {
				divergence[index] = field_x[index + 1];
			} else if (it->x == size_x - 1) {
				divergence[index] = - field_x[index];
			} else {
				divergence[index] = field_x[index + 1] - field_x[index];
			}
		} else {
			if (it->x == 0) {
				divergence[index] = field_x[index];
			} else if (it->x == size_x - 1) {
				divergence[index] = - field_x[index - 1];
			} else {
				divergence[index] = field_x[index] - field_x[index - 1];
			}
		}

		if (is_y_forward) {
			if (it->y == 0) {
				divergence[index] += field_y[index + size_x];
			} else if (it->y == size_y - 1) {
				divergence[index] += - field_y[index];
			} else {
				divergence[index] += field_y[index + size_x] - field_y[index];
			}
		} else {
			if (it->y == 0) {
				divergence[index] += field_y[index];
			} else if (it->y == size_y - 1){
				divergence[index] += - field_y[index - size_x];
			} else {
				divergence[index] += field_y[index] - field_y[index - size_x];
			}
		}
	}
}


/* Calculates anisotropic laplacian with homogeneous Neumann boundary conditions on the
 * image boundary (not the inpainting domain boundary), i.e. the gradient is zero.
 */
inline void PatchNonLocalPoisson::calculate_anisotropic_laplacian(double *image,
														   const double *coefficients,
														   const FixedMask &mask)
{
	int size_x = mask.get_size_x();
	int size_y = mask.get_size_y();

	// NOTE: The boundary is the border of the image (rows 0 and height-1 and columns 0 and width-1).

	vector<Point> points = mask.get_masked_points();
	double* buffer = new double[points.size()];

	double a_top, a_bottom, a_right, a_left;	// coefficients
	double u_top, u_bottom, u_right, u_left;	// color values

	// calculate anisotropic laplacian and store it in the buffer
	for (unsigned int i = 0; i < points.size(); i++) {
		Point p = points[i];
		int index = size_x * p.y + p.x;

		if (p.y != 0) {
			a_top = coefficients[index - size_x];
			u_top = image[index - size_x];
		} else {
			a_top = 0;
			u_top = 0;
		}

		if (p.x != 0) {
			a_left =  coefficients[index - 1];
			u_left = image[index - 1];
		} else {
			a_left = 0;
			u_left = 0;
		}

		if (p.y != size_y - 1) {
			a_bottom = coefficients[index];	// at the center
			u_bottom = image[index + size_x];
		} else {
			a_bottom = 0;
			u_bottom = 0;
		}

		if (p.x != size_x - 1) {
			a_right = coefficients[index];	// at the center
			u_right = image[index + 1];
		} else {
			a_right = 0;
			u_right = 0;
		}

		// the four combinations averaged
		buffer[i] = u_top * a_top + u_left * a_left + u_bottom * a_bottom + u_right * a_right -
						image[index] * (a_top + a_left + a_bottom + a_right);
	}

	// copy data from the buffer to the inpainting domain
	for (unsigned int i = 0; i < points.size(); i++) {
		Point p = points[i];
		int index = size_x * p.y + p.x;
		image[index] = buffer[i];
	}

	delete[] buffer;
}


/**
 * Stores image data channel by channel.
 */
inline double** PatchNonLocalPoisson::split_image_into_channels(const FixedImage<float> &image)
{
	int size_x = image.get_size_x();
	int size_y = image.get_size_y();
	int number_of_channels = image.get_number_of_channels();

	// allocate memory
	double** splitted_image = new double*[number_of_channels];
	for (int ch = 0; ch < number_of_channels; ch++) {
		splitted_image[ch] = new double[size_x * size_y];
	}

	const float* image_data = image.raw();

	// split image data
	for (int i = 0; i < size_x * size_y; i++) {
		for (int ch = 0; ch < number_of_channels; ch++) {
			splitted_image[ch][i] = image_data[i * number_of_channels + ch];
		}
	}

	return splitted_image;
}


inline void PatchNonLocalPoisson::conjugte_gradient(double *perturbation,
											 const double *a1,
											 const double *a2,
											 const double *b,
											 const FixedMask& mask)
{
	int size_x = mask.get_size_x();
	int size_y = mask.get_size_y();
	int pixels_amount = size_x * size_y;

	vector<Point> points = mask.get_masked_points();

	// allocate memory
	double *al_image = new double[pixels_amount];	// 'al_image' stands for the image with applied anisotropic laplacian
	double *d = new double[pixels_amount]();	// NOTE: default-initialized to 0.0
	double *r = new double[points.size()]();	// NOTE: default-initialized to 0.0

	/// initializations

	// r = b - A.x ; nr = |r| ; d = r
	memcpy(al_image, perturbation, pixels_amount * sizeof(double));
	calculate_anisotropic_laplacian(al_image, a2, mask);

	double nr0 = 0.0;
	for (unsigned int i = 0; i < points.size(); i++)	{
		Point p = points[i];
		int index = size_x * p.y + p.x;
		r[i] = b[index] - a1[index] * perturbation[index] + al_image[index];
		d[index] = r[i];
		nr0 += r[i] * r[i];
	}

	/// main loop
	int iteration = 0;
	double nr = nr0;
	bool stop = !(nr0 > 0);

	while (!stop) {
		// alpha = |r| / <r,q> ; q = A.d
		memcpy(al_image, d, pixels_amount * sizeof(double));
		calculate_anisotropic_laplacian(al_image, a2, mask);
		double rdotq = 0.0;

		for (unsigned int i = 0; i < points.size(); i++)	{
			Point p = points[i];
			int index = size_x * p.y + p.x;
			rdotq += r[i] * (a1[index] * d[index] - al_image[index]);
		}
		double alpha = nr/rdotq;

		// x = x + alpha*d
		for (unsigned int i = 0; i < points.size(); i++) {
			Point p = points[i];
			int index = size_x * p.y + p.x;
			perturbation[index] += alpha * d[index];
		}

		// r = b - A.u ; nr = |r|
		double nr_old = nr;
		memcpy(al_image, perturbation, pixels_amount * sizeof(double));
		calculate_anisotropic_laplacian(al_image, a2, mask);
		nr = 0.0;
		for (unsigned int i = 0; i < points.size(); i++) {
			Point p = points[i];
			int index = size_x * p.y + p.x;
			r[i] = b[index] - a1[index] * perturbation[index] + al_image[index];
			nr += r[i] * r[i];
		}

		// beta = nr/nr_old ; d = r + beta.d
		double beta = nr/nr_old;
		for (unsigned int i = 0; i < points.size(); i++) {
			Point p = points[i];
			int index = size_x * p.y + p.x;
			d[index] = r[i] + beta * d[index];
		}

		// check stopping
		stop = ((nr < _conjugate_gradient_tolerance*nr0) || (iteration >= _conjugate_gradient_iterations_limit));
		iteration++;
	}

#ifdef DBG_OUTPUT
//	printf("\t\t\tconj.grad. ended with residue %g after %d iterations\n", nr, iteration);
#endif

	// free memory
	delete[] al_image;
	delete[] d;
	delete[] r;
}
