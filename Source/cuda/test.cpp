#include "dmat.h"
#include <math.h>

void distanceCPU(const float*vg_a, int n_a, int k, float* d) {
	int i, j, att;

	for(i=0; i < n_a; i++) {
		for(j=0; j<n_a; j++) {

			if( i == j)
				d[i*n_a+j] = 0;
			else {
				float dis = 0;
				for(att = 0; att < k; att++)
					dis += (vg_a[i*k+att] - vg_a[j*k+att])*(vg_a[i*k+att] - vg_a[j*k+att]);
				dis = sqrt(dis);
				d[i*n_a+j] = dis;
			}
			
		}
	}
}

int main(int argc, const char* argv[]) {

	//void distance(	const float * vg_a, size_t pitch_a, size_t n_a,
	//			size_t k, float * d, size_t pitch_d );

	const float values[11][9] = { 	{0.181254, -0.585135, -0.60927, 0.191352, -0.402453, 0.14695, -0.282914, 0.519187, -0.00745311}, 
								{-0.114886, -0.286688, -0.00296578, -0.202359, -0.239975, 0.171379, -0.0611926, 0.45349, -0.0935648},
								{0.589659, -0.641492, 0.489091, 0.376823, -0.228953, 0.47934, 0.136344, -0.340885, 0.0764228},
								{-0.960732, 1, 0.293968, -0.0626467, -0.220027, -0.523963, -0.05336, 0.522304, 0.328645},
								{-0.84399, -0.543679, -0.362869, -0.253279, 0.225036, -0.789166, 0.253452, -0.685578, 0.311994}, 
								{-1, 0.359629, 1, -0.676075, -0.025972, 0.412234, 0.584265, -0.475419, -0.111974},
								{0.937092, 0.102344, -0.0282643, 1, 0.641631, -1, -0.192514, 0.601135, -1},
								{0.562969, 0.254497, -0.153238, 0.0401499, -0.649346, 0.358929, -1, 1, -0.179246}, 
								{0.544558, 0.71944, -0.69887, -0.864282, -0.503266, 0.189283, -0.21846, -0.768047, -0.0652458}, 
								{-0.334417, 0.159636, 0.154402, 0.0189855, 1, 0.1115, 0.465972, -0.662582, 0.942581},
								{0.438492, -0.538552, -0.0819845, 0.431333, 0.403324, 0.443514, 0.368407, -0.163604, -0.202159} };
	
	size_t k = 9;
	size_t pitch_a = k * sizeof(float);
	size_t n_a = 11;
	size_t pitch_d = n_a * sizeof(float);
	float* d = new float[n_a*n_a];


	float vg_a[99];

	for(int i=0; i < 11; i++)
		memcpy(vg_a+i*k, values[i], k*sizeof(float));

	for(int i = 0; i < 11; i++) {
		for(int j = 0; j < 9; j++) {
			//vg_a[i*9 + j] = values[i][j];
			printf("%10.5f", vg_a[i*9 + j]);
		}
		printf("\n");
	}

	distance(vg_a, pitch_a, n_a, k, 
				d, pitch_d);
	
	printf("GPU version:\n");
	for(int i = 0; i < n_a; i++) {
		for(int j = 0; j < n_a; j++) {
			printf("%10.5f", d[i*n_a + j]);
		}
		printf("\n");
	}

	printf("\nCPU version:\n");
	distanceCPU(vg_a, n_a, k, d);
	for(int i = 0; i < n_a; i++) {
		for(int j = 0; j < n_a; j++) {
			printf("%10.5f", d[i*n_a + j]);
		}
		printf("\n");
	}
	
	delete [] d; 

	return 0;
}