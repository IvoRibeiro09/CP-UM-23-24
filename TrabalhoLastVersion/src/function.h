#include "function.cpp"

void initialize();

double MeanSquaredVelocity();

double Kinetic();

double Potential();

void computeAccelerations();

double VelocityVerlet(double dt, int iter, FILE *fp);

void initializeVelocities();

double gaussdist();