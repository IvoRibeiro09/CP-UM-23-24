#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>
#include <cuda.h>
#include <chrono>

cudaEvent_t start, stop;

using namespace std;

void computeAccelerationsKernel(double *r, double *a, double *PE, int VECSIZE);

void computeAccelerations_plus_potential_GPU(double *r, double *a, double *PE, int VECSIZE);

  