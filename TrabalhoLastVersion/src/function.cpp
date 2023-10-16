// Number of particles
int N;

//  Lennard-Jones parameters in natural units!
double sigma = 1.;
double epsilon = 1.;
double m = 1.;
double kB = 1.;

double NA = 6.022140857e23;
double kBSI = 1.38064852e-23;  // m^2*kg/(s^2*K)

//  Size of box, which will be specified in natural units
double L;

//  Initial Temperature in Natural Units
double Tinit;  //2;
//  Vectors!
//
const int MAXPART=5001;
//  Position
double r[MAXPART][3];
//  Velocity
double v[MAXPART][3];
//  Acceleration
double a[MAXPART][3];
//  Force
double F[MAXPART][3];

// atom type
char atype[10];


//  Numerical recipes Gaussian distribution number generator
double gaussdist() {
    static bool available = false;
    static double gset;
    double fac, rsq, v1, v2;
    if (!available) {
        do {
            v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        
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
    
    for (i=0; i<N; i++) {
        /*
        for (j=0; j<3; j++) {
            //  Pull a number from a Gaussian Distribution
            v[i][j] = gaussdist();
            
        }*/
        v[i][0] = gaussdist();
        v[i][1] = gaussdist();
        v[i][2] = gaussdist();
    }
    
    // Vcm = sum_i^N  m*v_i/  sum_i^N  M
    // Compute center-of-mas velocity according to the formula above
    double vCM[3] = {0, 0, 0};
    
    for (i=0; i<N; i++) {
        /*
        for (j=0; j<3; j++) {
            
            vCM[j] += m*v[i][j];
            
        }*/
        vCM[0] += m*v[i][0];
        vCM[1] += m*v[i][1];
        vCM[2] += m*v[i][2];
    }
    
    /*
    for (i=0; i<3; i++) vCM[i] /= N*m;
    */
    vCM[0] /= N*m;
    vCM[1] /= N*m;
    vCM[2] /= N*m;


    //  Subtract out the center-of-mass velocity from the
    //  velocity of each particle... effectively set the
    //  center of mass velocity to zero so that the system does
    //  not drift in space!
    for (i=0; i<N; i++) {
        /*
        for (j=0; j<3; j++) {
            
            v[i][j] -= vCM[j];
            
        }*/
        v[i][0] -= vCM[0];
        v[i][1] -= vCM[1];
        v[i][2] -= vCM[2];
    }
    
    //  Now we want to scale the average velocity of the system
    //  by a factor which is consistent with our initial temperature, Tinit
    double vSqdSum, lambda;
    vSqdSum=0.;
    for (i=0; i<N; i++) {
        /*
        for (j=0; j<3; j++) {
            
            vSqdSum += v[i][j]*v[i][j];
            
        }*/
        vSqdSum += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
    }
    
    lambda = sqrt( 3*(N-1)*Tinit/vSqdSum);
    
    for (i=0; i<N; i++) {
        /*
        for (j=0; j<3; j++) {
            
            v[i][j] *= lambda;
            
        }
        */
        v[i][0] *= lambda;
        v[i][1] *= lambda;
        v[i][2] *= lambda;

    }
}

void initialize() {
    int n, p, i, j, k;
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
                    
                    r[p][0] = (i + 0.5)*pos;
                    r[p][1] = (j + 0.5)*pos;
                    r[p][2] = (k + 0.5)*pos;
                }
                p++;
            }
        }
    }
    
    // Call function to initialize velocities
    initializeVelocities();
    
    /***********************************************
     *   Uncomment if you want to see what the initial positions and velocities are
     printf("  Printing initial positions!\n");
     for (i=0; i<N; i++) {
     printf("  %6.3e  %6.3e  %6.3e\n",r[i][0],r[i][1],r[i][2]);
     }
     
     printf("  Printing initial velocities!\n");
     for (i=0; i<N; i++) {
     printf("  %6.3e  %6.3e  %6.3e\n",v[i][0],v[i][1],v[i][2]);
     }
     */
    
    
    
}   


//  Function to calculate the averaged velocity squared
double MeanSquaredVelocity() { 

    double velo = 0.0;
    
    for(int i=0; i<N ; i++){
        velo += v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2];
    }

    return velo = velo/N;

    /*
    double vx2 = 0;
    double vy2 = 0;
    double vz2 = 0;
    double v2;
    
    for (int i=0; i<N; i++) {
        
        vx2 = vx2 + v[i][0]*v[i][0];
        vy2 = vy2 + v[i][1]*v[i][1];
        vz2 = vz2 + v[i][2]*v[i][2];
        
    }
    v2 = (vx2+vy2+vz2)/N;
    
    
    //printf("  Average of x-component of velocity squared is %f\n",v2);
    return v2;
    */
}

//  Function to calculate the kinetic energy of the system
double Kinetic(){ //Write Function here!  

    double kin = 0.,v2;

    for(int i=0; i<N; i++){
        kin += m*(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2])/2.;
    }
    return kin;
    
    /*
    double v2, kin;
    
    kin =0.;
    for (int i=0; i<N; i++) {
        
        v2 = 0.;
        for (int j=0; j<3; j++) {
            
            v2 += v[i][j]*v[i][j];
            
        }
        kin += m*v2/2.;
        
    }
    
    //printf("  Total Kinetic Energy is %f\n",N*mvs*m/2.);
    return kin;
    */
}


// Function to calculate the potential energy of the system
double Potential() {
    double quot, r2, rnorm, term1, term2, Pot;
    int i, j, k;
    
    Pot=0.;
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            
            if (j!=i) {
                r2=0.;
                //r2 +=((r[i][0]-r[j][0])*(r[i][0]-r[j][0]) + (r[i][1]-r[j][1])*(r[i][1]-r[j][1]) + (r[i][2]-r[j][2])*(r[i][2]-r[j][2]));
                
                for (k=0; k<3; k++) {
                    r2 += (r[i][k]-r[j][k])*(r[i][k]-r[j][k]);
                }
                
                rnorm=sqrt(r2);
                quot=sigma/rnorm;
                double aux1 = quot * quot * quot * quot;
                double aux2 = quot * quot * quot;
                term1 = aux1 * aux1 * aux1;
                term2 = aux2 * aux2;
                
                Pot += 4*epsilon*(term1 - term2);
                
            }
        }
    }
    
    return Pot;
}



//   Uses the derivative of the Lennard-Jones potential to calculate
//   the forces on each atom.  Then uses a = F/m to calculate the
//   accelleration of each atom. 
void computeAccelerations() {
    int i, j, k;
    double f, rSqd;
    double rij[3]; // position of i relative to j
    
     // After the main loop, set all accelerations to zero
    for (k = 0; k < 3; k++) {
        for (i = 0; i < N; i++) {
            a[i][k] = 0.;
        }
    }
    
    for (i = 0; i < N-1; i++) {
        for (j = i+1; j < N; j++) {
            // initialize r^2 to zero
            rSqd = 0.;

            for (k = 0; k < 3; k++) {
                // component-by-component position of i relative to j
                rij[k] = r[i][k] - r[j][k];
                // sum of squares of the components
                rSqd += rij[k] * rij[k];
            }

            // From derivative of Lennard-Jones with sigma and epsilon set equal to 1 in natural units!
            //f = 24. * (2 * pow(rSqd, -7) - pow(rSqd, -4));
            double aux1 = rSqd * rSqd * rSqd;
            f = (48.-24. * aux1)/(aux1 * aux1 * rSqd);

            for (k = 0; k < 3; k++) {
                // from F = ma, where m = 1 in natural units!
                double aux = rij[k] * f;
                a[i][k] += aux;
                a[j][k] -= aux;
            }
        }
    }
}

// returns sum of dv/dt*m/A (aka Pressure) from elastic collisions with walls
double VelocityVerlet(double dt, int iter, FILE *fp) {
    int i, j, k;
    
    double psum = 0.;
    
    //  Compute accelerations from forces at current position
    // this call was removed (commented) for predagogical reasons
    //computeAccelerations();
    //  Update positions and velocity with current velocity and acceleration
    //printf("  Updated Positions!\n");
    for (i=0; i<N; i++) {
        for (j=0; j<3; j++) {
            r[i][j] += v[i][j]*dt + 0.5*a[i][j]*dt*dt;
            
            v[i][j] += 0.5*a[i][j]*dt;
        }
        //printf("  %i  %6.4e   %6.4e   %6.4e\n",i,r[i][0],r[i][1],r[i][2]);
    }
    //  Update accellerations from updated positions
    computeAccelerations();
    //  Update velocity with updated acceleration
    for (i=0; i<N; i++) {
        for (j=0; j<3; j++) {
            v[i][j] += 0.5*a[i][j]*dt;
        }
    }
    
    // Elastic walls
    for (i=0; i<N; i++) {
        for (j=0; j<3; j++) {
            if (r[i][j]<0.) {
                v[i][j] *=-1.; //- elastic walls
                psum += 2*m*fabs(v[i][j])/dt;  // contribution to pressure from "left" walls
            }
            if (r[i][j]>=L) {
                v[i][j]*=-1.;  //- elastic walls
                psum += 2*m*fabs(v[i][j])/dt;  // contribution to pressure from "right" walls
            }
        }
    }
    
    
    /* removed, uncomment to save atoms positions */
    /*for (i=0; i<N; i++) {
        fprintf(fp,"%s",atype);
        for (j=0; j<3; j++) {
            fprintf(fp,"  %12.10e ",r[i][j]);
        }
        fprintf(fp,"\n");
    }*/
    //fprintf(fp,"\n \n");
    
    return psum/(6*L*L);
}