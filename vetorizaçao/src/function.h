#include "function.cpp"

void initialize();

double MeanSquaredVelocity();

double Kinetic();

void Potential();

void computeAccelerations();

//void computeAccelerations_plus_potencial();

void clearAmatrix();

double VelocityVerlet(double dt, int iter, FILE *fp);

void initializeVelocities();

double gaussdist();
