#include "function.cpp"
#include <cmath>

void initialize();

void MeanSquaredVelocity();

void Kinetic();

void Potential();

void MeanSquaredVelocity_and_Kinetic();

void make_m_arr();

void computeAccelerations();

void computeAccelerations_plus_potential();

void clearAmatrix();

double VelocityVerlet(double dt, int iter, FILE *fp);

void initializeVelocities();

double gaussdist();
