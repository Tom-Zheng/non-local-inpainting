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

#include "a_image_updating.h"

AImageUpdating::AImageUpdating()
{
	_gaussian_sigma = 1.0;
	_patch_size = Shape(7, 7);
}


AImageUpdating::AImageUpdating(Shape patch_size, float gaussian_sigma)
{
	_gaussian_sigma = gaussian_sigma;
	_patch_size = patch_size;
}


AImageUpdating::~AImageUpdating()
{

}


float AImageUpdating::get_gaussian_sigma()
{
	return _gaussian_sigma;
}


void AImageUpdating::set_gaussian_sigma(float gaussian_sigma)
{
	if (_gaussian_sigma != gaussian_sigma) {
		_gaussian_sigma = gaussian_sigma;

		if (_patch_weighting.is_not_empty()) {
			_patch_weighting = Image<float>();
		}
	}
}


Shape AImageUpdating::get_patch_size()
{
	return _patch_size;
}


void AImageUpdating::set_patch_size(Shape patch_size)
{
	if (_patch_size != patch_size) {
		_patch_size = patch_size;

		if (_patch_weighting.is_not_empty()) {
			_patch_weighting = Image<float>();
		}
	}
}




