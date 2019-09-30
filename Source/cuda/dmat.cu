#include <stdio.h>
#include <cublas.h>

#include "dmat.h"

#define NUM_THREADS 32

// Space for the vector data
__constant__ float * distance_vg_a_d;

// Space for the resulting distance
__device__ float * distance_d_d;

void checkCudaError(const char * msg) {
	cudaError_t err = cudaGetLastError();
	if(cudaSuccess != err) {
		if(msg != NULL) {
			printf("[WARNING] %s\n", msg);
		}
		printf("[ERROR] %s : %s\n", msg, cudaGetErrorString(err));
		exit(1);
	}
}

__global__ void euclidean_kernel_same(	const float * vg_a, size_t pitch_a, size_t n_a, size_t k,
										float * d, size_t pitch_d )
{
	size_t x = blockIdx.x, y = blockIdx.y;

	if((x == y) && (x < n_a) && (threadIdx.x == 0))
		d[y * pitch_d + x] = 0.0;
  
	// If all element is to be computed
	if(y < n_a && x < y) {
		__shared__ float temp[NUM_THREADS];    

		temp[threadIdx.x] = 0.0;
    
		for(size_t offset = threadIdx.x; offset < k; offset += NUM_THREADS) {
			float t = abs(vg_a[x * pitch_a + offset] - vg_a[y * pitch_a + offset]);
			temp[threadIdx.x] += (t * t);
		}
    
		// Sync with other threads
		__syncthreads();
    
		// Reduce
		for(size_t stride = blockDim.x >> 1; stride > 0; stride >>= 1) {
			if(threadIdx.x < stride)
				temp[threadIdx.x] += temp[threadIdx.x + stride];
			__syncthreads();
		}
	    
		// Write to global memory
		if(threadIdx.x == 0) {
			float s = sqrt(temp[0]);
			d[y * pitch_d + x] = s;
			d[x * pitch_d + y] = s;
		}
	}
}

void distance_device(	const float * vg_a_d, size_t pitch_a, size_t n_a, size_t k,
						float * d_d, size_t pitch_d) {

	dim3 block(NUM_THREADS, 1, 1);
	dim3 grid(n_a, n_a, 1);

	size_t fbytes = sizeof(float);

	pitch_a /= fbytes;
	pitch_d /= fbytes;

	euclidean_kernel_same<<<grid, block>>>(	vg_a_d, pitch_a, n_a, k,
											d_d, pitch_d);
}

void distanceGPU(	const float *vg_a, size_t pitch_a, size_t n_a, size_t k, 
				float * d, size_t pitch_d ) {
	
	size_t pitch_a_d, pitch_d_d;
	
	// Allocate space for the vectors and distances on the gpu
	cudaMallocPitch((void**)&distance_vg_a_d, &pitch_a_d, k * sizeof(float), n_a);
	cudaMemcpy2D(distance_vg_a_d, pitch_a_d, vg_a, pitch_a, k * sizeof(float), n_a, cudaMemcpyHostToDevice);

	cudaMallocPitch((void**)&distance_d_d, &pitch_d_d, n_a * sizeof(float), n_a);
	
	checkCudaError("distance function : malloc and memcpy");
    
	distance_device(distance_vg_a_d, pitch_a_d, n_a, k,
					distance_d_d, pitch_d_d);
	
	checkCudaError("distance function : kernel invocation");

	// Copy the result back to cpu land now that gpu work is done
	cudaMemcpy2D(d, pitch_d, distance_d_d, pitch_d_d, n_a * sizeof(float), n_a, cudaMemcpyDeviceToHost);
	checkCudaError("distance function : memcpy");
    
	// Free allocated space
	cudaFree(distance_vg_a_d);
	cudaFree(distance_d_d);
}