#include <stdio.h>
#include <cstdlib>
#include <iostream>
#include <sys/time.h>
#include <cuda.h>
#include <chrono>
#include <stdlib.h>
#include <math.h>
#include <string.h>

cudaEvent_t start, stop;

using namespace std;

void checkCUDAError (const char *msg) {
	cudaError_t err = cudaGetLastError();
	if( cudaSuccess != err) {
		cerr << msg << "\n" << cudaGetErrorString( err) << endl;
		exit(-1);
	}
}

void startKernelTime (void) {
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start);
}

void stopKernelTime (void) {
	cudaEventRecord(stop);

	cudaEventSynchronize(stop);
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);

	cout << "\n" << "Basic profiling: " << milliseconds << " ms have elapsed for the kernel execution" << "\n" << endl;
}