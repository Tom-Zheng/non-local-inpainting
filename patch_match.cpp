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

#include "patch_match.h"

PatchMatch::PatchMatch()
{
	_iteration_count = 10;
	_search_window_size = -1;
	_random_shots_limit = 20;
	_distance_calculation = 0;
	_max_random_shots_count = 0;
}


PatchMatch::PatchMatch(APatchDistance *distance_calculation)
{
	_iteration_count = 10;
	_search_window_size = -1;
	_random_shots_limit = 20;
	_distance_calculation = distance_calculation;
	_max_random_shots_count = 0;
}


PatchMatch::PatchMatch(int iteration_count, int random_shots_limit, int search_window_size)
{
	_iteration_count = iteration_count;
	_search_window_size = search_window_size;
	_random_shots_limit = random_shots_limit;
	_distance_calculation = 0;
	_max_random_shots_count = 0;
}


PatchMatch::PatchMatch(APatchDistance *distance_calculation,
				int iteration_count,
				int random_shots_limit,
				int search_window_size)
{
	_iteration_count = iteration_count;
	_search_window_size = search_window_size;
	_random_shots_limit = random_shots_limit;
	_distance_calculation = distance_calculation;
	_max_random_shots_count = 0;
}

#ifdef _OPENMP

/**
 * Estimates NNF using the given initial nearest neighbors field.
 * Uses OpenMP for parallelization.
 *
 * @param initial_field Initial nearest neighbors field. Null pointer causes random initialization.
 */
Image<Point> PatchMatch::calculate(FixedImage<float> source,
									FixedMask source_mask,
									FixedImage<float> target,
									FixedMask target_mask,
									Image<Point> initial_field)
{
	if ((!initial_field.is_empty() && initial_field.get_size() != target.get_size()) ||
			(source.get_size() != source_mask.get_size()) ||
			(target.get_size() != target_mask.get_size()) ||
			!_distance_calculation) {
		return Image<Point>();
	}

	// Initialize distance calculation
	_distance_calculation->initialize(source, target);

	Shape target_shape = target.get_size();
	Shape source_shape = source.get_size();

	// Allocate memory for nearest neighbors and distances
	// NOTE: we need two buffers for both distances and neighbors to avoid data access conflicts for adjacent threads.
	//		 Threads with odd indices work with *_odd buffers, while threads with even indices work with *_even.
	Image<float> distances_odd(target_shape.size_x, target_shape.size_y, numeric_limits<float>::max());
	Image<float> distances_even(target_shape.size_x, target_shape.size_y, numeric_limits<float>::max());
	Image<Point> neighbors_odd(target_shape.size_x, target_shape.size_y, Point(-1, -1));
	Image<Point> neighbors_even(target_shape.size_x, target_shape.size_y, Point(-1, -1));

	// Build masked points cache for speedup
	vector<Point> target_points = target_mask.get_masked_points();

	// Use given nearest neighbor field (NNF) or initialize NNF at random.
	if (!initial_field.is_empty() && initial_field.get_size() == target_shape) {
		// Reinitialize shifts pointing outside the target region and calculate distances
		for (uint i = 0; i < target_points.size(); i++) {
			Point p = target_points[i];
			Point neighbor = initial_field(p);

			int number_of_tries = 0;
			while (!source_mask.test(neighbor.x, neighbor.y) && number_of_tries < _random_shots_limit) {
				neighbor.x = rand() % source_shape.size_x;
				neighbor.y = rand() % source_shape.size_y;
				number_of_tries++;
			}

			if (source_mask.test(neighbor.x, neighbor.y)) {
				neighbors_odd(p) = neighbor;
				neighbors_even(p) = neighbor;

				float distance = _distance_calculation->calculate(neighbor, p);
				distances_odd(p) = distance;
				distances_even(p) = distance;
			}
		}
	} else {
		// Initialize shifts at random and calculate distances
		for (uint i = 0; i < target_points.size(); i++) {
			Point p = target_points[i];
			Point neighbor = Point(-1, -1);

			int number_of_tries = 0;
			while ((neighbor.x < 0 || !source_mask.get(neighbor.x, neighbor.y)) && number_of_tries < _random_shots_limit) {
				neighbor.x = rand() % source_shape.size_x;
				neighbor.y = rand() % source_shape.size_y;
				number_of_tries++;
			}

			if (source_mask.test(neighbor.x, neighbor.y)) {
				neighbors_odd(p) = neighbor;
				neighbors_even(p) = neighbor;

				float distance = _distance_calculation->calculate(neighbor, p);
				distances_odd(p) = distance;
				distances_even(p) = distance;
			}
		}
	}

	int inpainting_domain_width = target_mask.bounding_box_bottom_right().x - target_mask.bounding_box_top_left().x + 1;

	// Base seed for random number generator
	uint seed = time(NULL);

	// NOTE: each thread should get the number of target points not less then doubled inpainting domain width.
	//       In this case we can safely copy data from one buffer to another after each iteration.
	#pragma omp parallel firstprivate(seed) num_threads( min(omp_get_max_threads(), (int)target_points.size() / (int)(2 * inpainting_domain_width)) )
	{	// === start of parallel block ===

		// Get thread-specific data
#ifdef _OPENMP
		int thread_id = omp_get_thread_num();
		int number_of_threads = omp_get_num_threads();
#else 
		int thread_id = 0;
		int number_of_threads = 1;
#endif

		// Specify seed for each thread
		seed += thread_id;

		int chunk_size = target_points.size() / number_of_threads;

		// Initialize appropriate shortcuts for buffers
		Image<Point> *my_neighbors;
		Image<Point> *other_neighbors;
		Image<float> *my_distances;
		Image<float> *other_distances;
		if ( thread_id % 2 != 0 ) {
			my_neighbors = &neighbors_odd;
			other_neighbors = &neighbors_even;
			my_distances = &distances_odd;
			other_distances = &distances_even;
		} else {
			my_neighbors = &neighbors_even;
			other_neighbors = &neighbors_odd;
			my_distances = &distances_even;
			other_distances = &distances_odd;
		}

		// In each iteration, improve the NNF, by looping in scanline or reverse-scanline order.
		for (int iter = 0; iter < _iteration_count; iter++) {
			// Iterate forward in even iteration and backward in odd ones (indices depend on the thread id)
			int index_begin, index_end, shift;
			if ( iter % 2 == 0 ) {
				index_begin = chunk_size * thread_id;
				index_end = (thread_id < number_of_threads - 1) ? chunk_size * (thread_id + 1) : target_points.size();
				shift = -1;
			} else {
				index_begin = (thread_id < number_of_threads - 1) ? chunk_size * (thread_id + 1) - 1 : target_points.size() - 1;
				index_end = chunk_size * thread_id - 1;
				shift = 1;
			}

			for (int index = index_begin; index != index_end; index -= shift) {
				int x = target_points[index].x;
				int y = target_points[index].y;

				float distance = (*my_distances)(x, y);
				float original_distance = distance;
				Point neighbor(-1, -1);

				/// Propagation: Improve current guess by trying instead correspondences from left and above (below and right on odd iterations).
				if (target_mask.test(x + shift, y)) {
					Point candidate = (*my_neighbors)(x + shift, y);
					candidate.x -= shift;

					if (source_mask.test(candidate.x, candidate.y)) {
						// Check for improvement
						float candidate_distance = _distance_calculation->calculate(candidate, Point(x, y));
						if (candidate_distance < distance) {
							distance = candidate_distance;
							neighbor = candidate;
						}
					}
				}

				if (target_mask.test(x, y + shift)) {
					Point candidate = (*my_neighbors)(x, y + shift);
					candidate.y -= shift;

					if (source_mask.test(candidate.x, candidate.y)) {
						// Check for improvement
						float candidate_distance = _distance_calculation->calculate(candidate, Point(x, y));
						if (candidate_distance < distance) {
							distance = candidate_distance;
							neighbor = candidate;
						}
					}
				}

				if (neighbor.x < 0) {
					neighbor = (*my_neighbors)(x, y);
				}

				/// Random search: Improve current guess by searching in boxes of exponentially decreasing size around the current best guess.
				int max_window_size = (_search_window_size != -1) ? _search_window_size :
																	std::max(source_shape.size_x, source_shape.size_y);

				Point search_center = neighbor;
				for (int window_size = max_window_size; window_size >= 1; window_size /= 2) {
					// Limit sampling window
					int x_min = max(search_center.x - window_size, 0);
					int y_min = max(search_center.y - window_size, 0);
					int x_max = min(search_center.x + window_size + 1, (int)source_shape.size_x);
					int y_max = min(search_center.y + window_size + 1, (int)source_shape.size_y);

					// Sample
					Point candidate;
					for (int k = 0; k < _random_shots_limit; k++)
					{
						candidate.x = x_min + rand_r(&seed) % (x_max - x_min);
						candidate.y = y_min + rand_r(&seed) % (y_max - y_min);

						if (source_mask.test(candidate.x, candidate.y)) {
							// Check for improvement
							float candidate_distance = _distance_calculation->calculate(candidate, Point(x, y));
							if (candidate_distance < distance) {
								distance = candidate_distance;
								neighbor = candidate;
							}

							break;
						}
					}	// for (int k = 0; k < _random_shots_limit; k++)
				}	// for (int window_size = max_window_size; window_size >= 1; window_size /= 2)

				if (original_distance > distance) {
					(*my_distances)(x, y) = distance;
					(*my_neighbors)(x, y) = neighbor;
				}

			}	// for (ind = ind_begin; ind != ind_end; ind -= shift)


			#pragma omp barrier

			// Copy values at the front boundary of the chunk to the second buffer to allow information propagation to the next thread.
			// NOTE: we do not calculate the precise number of points that have to be copied, instead we copy at most N points,
			//       where N is the width of the inpainting domain's bounding box. In this way we can be sure that we copy everything that is needed (and maybe a bit more).
			if (iter < _iteration_count - 1) {
				int count = 0;
				for (int index = index_end + shift; (index != index_begin + shift) && (count < inpainting_domain_width); index += shift, count++) {
					Point p = target_points[index];
					(*other_distances)(p) = (*my_distances)(p);
					(*other_neighbors)(p) = (*my_neighbors)(p);
				}
			} else {
				// Synchronize buffers (only neighbors) in the end of the last iteration.
				for (int ind = index_end + shift; ind != index_begin + shift; ind += shift) {
					(*other_neighbors)(target_points[ind]) = (*my_neighbors)(target_points[ind]);
				}
			}

			#pragma omp barrier
		} // for (int i = 0; i < _iteration_count; i++) {
	} // === end of parallel block ===

	return neighbors_odd;
}

#else	// undefined _OPENMP

/**
 * Estimates NNF using the given initial nearest neighbors field.
 *
 * @param initial_field Initial nearest neighbors field. Null pointer causes random initialization.
 */
Image<Point> PatchMatch::calculate(FixedImage<float> source,
									FixedMask source_mask,
									FixedImage<float> target,
									FixedMask target_mask,
									Image<Point> initial_field)
{
	if ((!initial_field.is_empty() && initial_field.get_size() != target.get_size()) ||
			(source.get_size() != source_mask.get_size()) ||
			(target.get_size() != target_mask.get_size()) ||
			!_distance_calculation) {
		return Image<Point>();
	}

#ifdef METRICS
	drop_metrics();
#endif

	// Initialize distance calculation
	_distance_calculation->initialize(source, target);

	// Allocate memory for nearest neighbors and distances
	Shape target_shape = target.get_size();
	Shape source_shape = source.get_size();
	Image<float> distances(target_shape.size_x, target_shape.size_y, numeric_limits<float>::max());
	Image<Point> neighbors(target_shape.size_x, target_shape.size_y, Point(-1, -1));

	// Build masked points cache for speedup
	vector<Point> target_points = target_mask.get_masked_points();

	// Use given nearest neighbor field (NNF) or initialize NNF at random.
	if (!initial_field.is_empty() && initial_field.get_size() == target_shape) {
		// Reinitialize shifts pointing outside the target region and calculate distances
		for (uint i = 0; i < target_points.size(); i++) {
			Point p = target_points[i];
			Point neighbor = initial_field(p);

			int number_of_tries = 0;
			while (!source_mask.test(neighbor.x, neighbor.y) && number_of_tries < _random_shots_limit) {
				neighbor.x = rand() % source_shape.size_x;
				neighbor.y = rand() % source_shape.size_y;
				number_of_tries++;
			}

			if (source_mask.test(neighbor.x, neighbor.y)) {
				neighbors(p) = neighbor;

				float distance = _distance_calculation->calculate(neighbor, p);
				distances(p) = distance;
			}
		}
	} else {
		// Initialize shifts at random and calculate distances
		for (uint i = 0; i < target_points.size(); i++) {
			Point p = target_points[i];
			Point neighbor = Point(-1, -1);

			int number_of_tries = 0;
			while ((neighbor.x < 0 || !source_mask.get(neighbor.x, neighbor.y)) && number_of_tries < _random_shots_limit) {
				neighbor.x = rand() % source_shape.size_x;
				neighbor.y = rand() % source_shape.size_y;
				number_of_tries++;
			}

			if (source_mask.test(neighbor.x, neighbor.y)) {
				neighbors(p) = neighbor;

				float distance = _distance_calculation->calculate(neighbor, p);
				distances(p) = distance;
			}
		}
	}

#ifdef METRICS
	// metrics
	int metric_max_random_shots_count;
	int metric_propagations_count;
	double metric_total_distance;
#endif

	// In each iteration, improve the NNF, by looping in scanline or reverse-scanline order.
	for (int iter = 0; iter < _iteration_count; iter++) {

#ifdef METRICS
		// Clear metrics
		metric_propagations_count = 0;
		metric_max_random_shots_count = 0;
		metric_total_distance = 0.0;
#endif

		// Iterate forward in even iteration and backward in odd ones.
		int index, index_end, shift;
		if ( iter % 2 == 0 ) {
			index = 0;
			index_end = target_points.size();
			shift = -1;
		} else {
			index = target_points.size() - 1;
			index_end = -1;
			shift = 1;
		}

		for (; index != index_end; index -= shift) {
			int x = target_points[index].x;
			int y = target_points[index].y;

			float distance = distances(x, y);
			float original_distance = distance;
			Point neighbor(-1, -1);

			/// Propagation: Improve current guess by trying instead correspondences from left and above (below and right on odd iterations).
			if (target_mask.test(x + shift, y)) {
				Point candidate = neighbors(x + shift, y);
				candidate.x -= shift;

				if (source_mask.test(candidate.x, candidate.y)) {
					// Check for improvement
					float candidate_distance = _distance_calculation->calculate(candidate, Point(x, y));
					if (candidate_distance < distance) {
						distance = candidate_distance;
						neighbor = candidate;
					}
				}
			}

			if (target_mask.test(x, y + shift)) {
				Point candidate = neighbors(x, y + shift);
				candidate.y -= shift;

				if (source_mask.test(candidate.x, candidate.y)) {
					// Check for improvement
					float candidate_distance = _distance_calculation->calculate(candidate, Point(x, y));
					if (candidate_distance < distance) {
						distance = candidate_distance;
						neighbor = candidate;
					}
				}
			}

			if (neighbor.x < 0) {
				neighbor = neighbors(x, y);
			}
#ifdef METRICS
			else {
				metric_propagations_count++;
			}
#endif

			/// Random search: Improve current guess by searching in boxes of exponentially decreasing size around the current best guess.
			int max_window_size = (_search_window_size != -1) ? _search_window_size :
																std::max(source_shape.size_x, source_shape.size_y);

			Point search_center = neighbor;
			for (int window_size = max_window_size; window_size >= 1; window_size /= 2) {
				// Limit sampling window
				int x_min = max(search_center.x - window_size, 0);
				int y_min = max(search_center.y - window_size, 0);
				int x_max = min(search_center.x + window_size + 1, (int)source_shape.size_x);
				int y_max = min(search_center.y + window_size + 1, (int)source_shape.size_y);

				// Sample
				Point candidate;
				for (int k = 0; k < _random_shots_limit; k++)
				{
					candidate.x = x_min + rand() % (x_max - x_min);
					candidate.y = y_min + rand() % (y_max - y_min);

					if (source_mask.test(candidate.x, candidate.y)) {
						// Check for improvement
						float candidate_distance = _distance_calculation->calculate(candidate, Point(x, y));
						if (candidate_distance < distance) {
							distance = candidate_distance;
							neighbor = candidate;
						}

#ifdef METRICS
						if (k > metric_max_random_shots_count) {
							metric_max_random_shots_count = k;
						}
#endif

						break;
					}
				}	// for (int k = 0; k < _random_shots_limit; k++)
			}	// for (int window_size = max_window_size; window_size >= 1; window_size /= 2)

			if (original_distance > distance) {
				distances(x, y) = distance;
				neighbors(x, y) = neighbor;
			}

#ifdef METRICS
			metric_total_distance += min(original_distance, distance);
#endif
		}	// for (; ind != ind_end; ind -= shift)

#ifdef METRICS
		push_metrics(metric_propagations_count, metric_max_random_shots_count, metric_total_distance);
#endif
	}	// for (int iter = 0; iter < _iteration_count; iter++)

	return neighbors;
}

#endif	// #ifdef _OPENMP


/* getters, setters */
int PatchMatch::get_iteration_count()
{
	return _iteration_count;
}

void PatchMatch::set_iteration_count(int iteration_count)
{
	_iteration_count = iteration_count;
}

int PatchMatch::get_search_window_size()
{
	return _search_window_size;
}

void PatchMatch::set_search_window_size(int search_window_size)
{
	_search_window_size = search_window_size;
}

int PatchMatch::get_random_shots_limit()
{
	return _random_shots_limit;
}

void PatchMatch::set_random_shots_limit(int random_shots_limit)
{
	_random_shots_limit = random_shots_limit;
}

void PatchMatch::set_distance_calculation(APatchDistance *distance_calculation)
{
	_distance_calculation = distance_calculation;
}


int PatchMatch::get_max_random_shots_metric()
{
	return _max_random_shots_count;
}


vector<int> PatchMatch::get_propagations_per_iteration_metric()
{
	return _propagations_per_iteration;
}


vector<double> PatchMatch::get_total_distance_per_iteration()
{
	return _total_distance_per_iteration;
}


/* Private */

#ifdef METRICS
void PatchMatch::push_metrics(int propagations_count, int max_random_shots_count, double total_distance)
{
	// store metrics
	_propagations_per_iteration.push_back(propagations_count);
	_total_distance_per_iteration.push_back(total_distance);
	_max_random_shots_count = max(_max_random_shots_count, max_random_shots_count);
}


void PatchMatch::drop_metrics()
{
	// drop stored metrics
	_propagations_per_iteration.clear();
	_total_distance_per_iteration.clear();
	_max_random_shots_count = 0;
}
#endif
