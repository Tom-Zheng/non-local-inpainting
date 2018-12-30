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

#ifndef PATCH_MATCH_H_
#define PATCH_MATCH_H_

/// NOTE: if METRICS is not defined, metrics will not be collected (but still can be requested, even though
///		 the value will always remain the same).
//#define METRICS

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <limits>
#include "image.h"
#include "mask.h"
#include "point.h"
#include "a_patch_distance.h"
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

/**
 * Implements Patch-Match algorithm (see "PatchMatch: A Randomized
 * Correspondence Algorithm for Structural Image Editing" by Barnes et al.)
 */
class PatchMatch
{
public:
	PatchMatch();
	PatchMatch(APatchDistance *distance_calculation);
	PatchMatch(int iteration_count, int random_shots_limit = 20, int search_window_size = -1);
	PatchMatch(APatchDistance *distance_calculation,
				int iteration_count,
				int random_shots_limit = 20,
				int search_window_size = -1);

	// Estimates NNF using the given initial nearest neighbors field.
	Image<Point> calculate(FixedImage<float> source,
						   FixedMask source_mask,
						   FixedImage<float> target,
						   FixedMask target_mask,
						   Image<Point> initial_field = Image<Point>());

	/// getters and setters for parameters
	int get_iteration_count();
	void set_iteration_count(int iteration_count);
	int get_search_window_size();
	void set_search_window_size(int search_window_size);
	int get_random_shots_limit();
	void set_random_shots_limit(int random_shots_limit);
	void set_distance_calculation(APatchDistance *distance_calculation);

	/// metrics
	int get_max_random_shots_metric();
	vector<int> get_propagations_per_iteration_metric();
	vector<double> get_total_distance_per_iteration();

private:
	APatchDistance *_distance_calculation;
	int _iteration_count;
	int _search_window_size;
	int _random_shots_limit;
	// metrics
	vector<int> _propagations_per_iteration;
	vector<double> _total_distance_per_iteration;
	int _max_random_shots_count;

#ifdef METRICS
	void push_metrics(int propagations_count, int max_random_shots_count, double total_distance);
	void drop_metrics();
#endif
};

#endif /* PATCH_MATCH_H_ */
