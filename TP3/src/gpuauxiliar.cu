#include "gpuauxiliar.h"
#define NUM_BLOCKS 128
#define NUM_THREADS_PER_BLOCK 256
#define SIZE NUM_BLOCKS*NUM_THREADS_PER_BLOCK

using namespace std;

__global__ void computeAccelerationsKernel(double *r, double *a, double *PE, int VECSIZE) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < VECSIZE) {
        double a0 = 0.0, a1 = 0.0, a2 = 0.0;
        double ri0 = r[3 * idx], ri1 = r[3 * idx + 1], ri2 = r[3 * idx + 2];

        for (int j = 0; j < VECSIZE; j++) {
            if (j != idx) {
                double M0 = ri0 - r[3 * j], M1 = ri1 - r[3 * j + 1], M2 = ri2 - r[3 * j + 2];
                double rSqd = M0 * M0 + M1 * M1 + M2 * M2;

                double aux = rSqd * rSqd * rSqd;
                double term2 = 1. / aux;
                double f = (48. - 24. * aux) / (aux * aux * rSqd);
                *PE += term2 * (term2 - 1.);

                double aux0 = M0 * f;
                double aux1 = M1 * f;
                double aux2 = M2 * f;

                a0 += aux0;
                a1 += aux1;
                a2 += aux2;
            }
        }

        atomicAdd(&a[3 * idx], a0);
        atomicAdd(&a[3 * idx + 1], a1);
        atomicAdd(&a[3 * idx + 2], a2);
    }
}

void computeAccelerations_plus_potential_GPU(double *r, double *a, double *PE, int VECSIZE) {
    *PE = 0.0;

    double *d_r, *d_a, *d_PE;

    cudaMalloc((void **)&d_r, sizeof(double) * 3 * VECSIZE);
    cudaMalloc((void **)&d_a, sizeof(double) * 3 * VECSIZE);
    cudaMalloc((void **)&d_PE, sizeof(double));

    cudaMemcpy(d_r, r, sizeof(double) * 3 * VECSIZE, cudaMemcpyHostToDevice);
    cudaMemcpy(d_a, a, sizeof(double) * 3 * VECSIZE, cudaMemcpyHostToDevice);
    cudaMemcpy(d_PE, PE, sizeof(double), cudaMemcpyHostToDevice);

    int blockSize = 256;
    int gridSize = (VECSIZE + blockSize - 1) / blockSize;

    computeAccelerationsKernel<<<gridSize, blockSize>>>(d_r, d_a, d_PE, VECSIZE);

    cudaMemcpy(a, d_a, sizeof(double) * 3 * VECSIZE, cudaMemcpyDeviceToHost);
    cudaMemcpy(PE, d_PE, sizeof(double), cudaMemcpyDeviceToHost);

    cudaFree(d_r);
    cudaFree(d_a);
    cudaFree(d_PE);
}