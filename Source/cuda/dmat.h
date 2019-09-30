#ifndef DMAT_H__
#define DMAT_H__

#include<stdio.h>
#include<string.h>

/*
Calculate the EUCLIDEAN distance matrix for vectors in group a and vectors in group b.
To calculate NORMALIZED EUCLIDEAN distance, normalise the input matrix first.
Ref: R gputools package

The format for vg_a is as follows:
There are n_* vectors, each of dimensionality k.  
They are stored in row major order with a row (or pitch) being pitch_* bytes.
The calculated distances are stored in d such that the distance between
vectors indexed a and b is located in d[b * pitch_d / sizeof(float) + a].
The user is responsible for allocating storage for d.  
It should be at least n_a * n_b * sizeof(float) bytes.  
The pitch_d argument is the same as for vg_*.
This function may run on the CPU or GPU depending on the size and
number of vectors.  Good data alignment will increase performance.
*/

void distanceGPU(	const float *vg_a, size_t pitch_a, size_t n_a, size_t k,
				float * d, size_t pitch_d );

void distance_device(	const float * vg_a_d, size_t pitch_a, size_t n_a, size_t k,
						float * d_d, size_t pitch_d );

#endif