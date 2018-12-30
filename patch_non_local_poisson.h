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

#ifndef PATCH_NON_LOCAL_POISSON_H_
#define PATCH_NON_LOCAL_POISSON_H_

#include <cstring>
#include "a_image_updating.h"
#include "gradient.h"

using namespace std;

/**
 * Implements Non-Local Poisson image updating scheme.
 */
class PatchNonLocalPoisson : public AImageUpdating
{
public:
	PatchNonLocalPoisson();
	PatchNonLocalPoisson(Shape patch_size,
						 float gaussian_sigma,
						 float lambda,
						 float conjugate_gradient_tolerance = 1e-10,
						 int conjugate_gradient_iterations_limit = 1000);

	virtual ~PatchNonLocalPoisson() {}

	virtual double update(Image<float> image,
						  Image<float> original_image,
						  FixedMask inpainting_domain,
						  FixedMask extended_inpainting_domain,
						  FixedImage<Point> nnf,
						   FixedImage<float> confidence_mask);

private:
	float _lambda;
	float _conjugate_gradient_tolerance;
	int _conjugate_gradient_iterations_limit;

	inline void calculate_pde_coefficients(const FixedImage<float> &image,
										   const FixedImage<float> &gradient,
										   const FixedMask &inpainting_domain,
										   const FixedImage<Point> &nnf,
										   const FixedImage<float> &confidence_mask,
										   double *a1,
										   double *a2,
										   double **f1,
										   double **f2);

	// Calculate the divergence with Neumann boundary conditions.
	inline void calculate_divergence(const double *field_x,
								     const double *field_y,
								     const FixedMask &mask,
								     bool is_x_forward,
								     bool is_y_forward,
								     double *out_divergence);

	inline void calculate_anisotropic_laplacian(double *image,
										 const double *coefficients,
										 const FixedMask &mask);

	inline double** split_image_into_channels(const FixedImage<float> &image);

	inline void conjugte_gradient(double *perturbation,
								  const double *a1,
								  const double *a2,
								  const double *b,
								  const FixedMask& mask);
};




#endif /* PATCH_NON_LOCAL_POISSON_H_ */
