//  Lennard-Jones parameters in natural units!
const int N = 2160;// Number of particles
const int ep_8 = 8;
const int MAXPART = 15000;
const int limit = 6480;
const int limitM1 = 6477;

double NA = 6.022140857e23;
double kBSI = 1.38064852e-23;  // m^2*kg/(s^2*K)
double m = 1.;
double kB = 1.;

double PE;
double mvs;
double KE;

//  Size of box, which will be specified in natural units
double L;

//  Initial Temperature in Natural Units
double Tinit;  //2;

//  Vectors!
//  Position
double r[MAXPART];
//  Velocity
double v[MAXPART];
//  Acceleration
double a[MAXPART];
//  Force
double F[MAXPART];
// atom type
char atype[10];

//  Numerical recipes Gaussian distribution number generator
double gaussdist() {
    static bool available = false;
    static double gset;
    double fac, rsq, v1, v2;
    if (!available) {
        do {
            v1 = 2. * rand() / double(RAND_MAX) - 1.;
            v2 = 2. * rand() / double(RAND_MAX) - 1.;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1. || rsq == 0.);
        
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        available = true;
        
        return v2*fac;
    } else {
        
        available = false;
        return gset;
        
    }
}

void initializeVelocities() {
    int i;
    for (i=0; i < limit;) {
        v[i++] = gaussdist();
        v[i++] = gaussdist();
        v[i++] = gaussdist();
    }
    
    // Vcm = sum_i^N  m*v_i/  sum_i^N  M
    // Compute center-of-mas velocity according to the formula above
    double vCM[3] = {0, 0, 0};
    
    for (i=0; i < limit;) {
        vCM[0] += m*v[i++];
        vCM[1] += m*v[i++];
        vCM[2] += m*v[i++];
    }
    
    vCM[0] /= N*m;
    vCM[1] /= N*m;
    vCM[2] /= N*m;
    //  Subtract out the center-of-mass velocity from the
    //  velocity of each particle... effectively set the
    //  center of mass velocity to zero so that the system does
    //  not drift in space!
    for (i=0; i < limit;) {
        v[i++] -= vCM[0];
        v[i++] -= vCM[1];
        v[i++] -= vCM[2];
    }
    
    //  Now we want to scale the average velocity of the system
    //  by a factor which is consistent with our initial temperature, Tinit
    double vSqdSum, lambda;
    vSqdSum=0.;
    for (i=0; i < limit;i += 3) {
        vSqdSum += (v[i]*v[i] + v[i + 1]*v[i + 1] + v[i + 2]*v[i + 2]);
    }
    
    lambda = sqrt( 3*(N-1)*Tinit/vSqdSum);
    
    for (i=0; i < limit;) {
        v[i++] *= lambda;
        v[i++] *= lambda;
        v[i++] *= lambda;
    }
}

void initialize() {
    int i, j, n, p, k;
    double pos;
    
    // Number of atoms in each direction
    n = int(ceil(pow(N, 1.0/3)));
    
    //  spacing between atoms along a given direction
    pos = L / n;
    
    //  index for number of particles assigned positions
    p = 0;
    //  initialize positions
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            for (k=0; k<n; k++) {
                if (p<N) {
                    
                    r[p*3+0] = (i + 0.5)*pos;
                    r[p*3+1] = (j + 0.5)*pos;
                    r[p*3+2] = (k + 0.5)*pos;
                }
                p++;
            }
        }
    }
    
    // Call function to initialize velocities
    initializeVelocities();

}   


//  Function to calculate the averaged velocity squared
//  Function to calculate the kinetic energy of the system
void MeanSquaredVelocity_and_Kinetic(){
    int i;
    double velo = 0., v2, kin;
    
    for(i=0; i < limit; i += 3){
        v2 = v[i]*v[i] + v[i + 1]*v[i + 1] + v[i + 2]*v[i + 2];
        velo += v2;
        kin += m*v2/2.;
    }
    KE = kin;
    mvs = velo/N;
}

// Function to calculate the potential energy of the system
//   Uses the derivative of the Lennard-Jones potential to calculate
//   the forces on each atom.  Then uses a = F/m to calculate the
//   accelleration of each atom. 
void computeAccelerations_plus_potential(){
    int i, j;
    double f, rSqd, term2, ri0, ri1, ri2, M0, M1, M2, aux, aux0, aux1, aux2, a0, a1, a2;
    PE = 0.;
    for (i = 0; i < limitM1; i += 3) { // loop over all distinct pairs i, j
        a0 = 0.0, a1 = 0.0, a2 = 0.0;

        ri0 = r[i];
        ri1 = r[i + 1];
        ri2 = r[i + 2];

        for (j = i + 3; j < limit; j += 3) {
            M0 = ri0 - r[j], M1 = ri1 - r[j + 1], M2 = ri2 - r[j + 2];

            rSqd = M0 * M0 + M1 * M1 + M2 * M2;
        
            aux = rSqd * rSqd * rSqd;
            term2 = 1. / aux;
            f = (48. - 24. * aux) / (aux * aux * rSqd);
            PE += ep_8 * term2 * (term2 - 1.);  

            aux0 = M0 * f;
            aux1 = M1 * f;
            aux2 = M2 * f;
            
            a0 += aux0;
            a1 += aux1;
            a2 += aux2;
            
            a[j] -= aux0, a[j + 1] -= aux1, a[j + 2] -= aux2;
        }
        a[i] += a0;
        a[i + 1] += a1;
        a[i + 2] += a2;
    }
}

//set all positions in a array to zero
void clearAmatrix(){
    int i;
    for (i = 0; i < limit;) {  // set all accelerations to zero
        a[i++] = 0.;
        a[i++] = 0.;
        a[i++] = 0.;

        a[i++] = 0.;
        a[i++] = 0.;
        a[i++] = 0.;
    }
}

// returns sum of dv/dt*m/A (aka Pressure) from elastic collisions with walls
double VelocityVerlet(double dt, int iter, FILE *fp) {
    double psum = 0.;
    int i, j;
    
    //  Compute accelerations from forces at current position
    // this call was removed (commented) for predagogical reasons
    //computeAccelerations();
    //  Update positions and velocity with current velocity and acceleration
    //printf("  Updated Positions!\n");
    for (i=0; i < limit; i +=3) {
        for (j=0; j<3; j++) {
            r[i+j] += v[i+j]*dt + 0.5*a[i+j]*dt*dt;
            
            v[i+j] += 0.5*a[i+j]*dt;
        }
        //printf("  %i  %6.4e   %6.4e   %6.4e\n",i,r[i][0],r[i][1],r[i][2]);
    }
    //  Update accellerations from updated positions
    clearAmatrix();
    computeAccelerations_plus_potential();
    //  Update velocity with updated acceleration
    for (i=0; i < limit; i +=3) {
        for (j=0; j<3; j++) {
            v[i+j] += 0.5*a[i+j]*dt;
            // Elastic walls
            if (r[i+j]<0.) {
                v[i+j] *=-1.; //- elastic walls
                psum += 2*m*fabs(v[i+j])/dt;  // contribution to pressure from "left" walls
            }
            if (r[i+j]>=L) {
                v[i+j]*=-1.;  //- elastic walls
                psum += 2*m*fabs(v[i+j])/dt;  // contribution to pressure from "right" walls
            }
        }
    }
    
    return psum/(6*L*L);
}

