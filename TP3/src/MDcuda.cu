#include "MDcuda.h"

#define NUM_THREADS_PER_BLOCK 12

// Number of particles
const int N = 5000;
// Vector's SIZE and vectorś size minus one
const int VECSIZE = 3 * N;
//const int VECSIZEM1 = 3 * (N-1);

__device__ int NCuda;


double PE, KE, mvs;
double NA = 6.022140857e23;
double kBSI = 1.38064852e-23;  // m^2*kg/(s^2*K)

//  Size of box, which will be specified in natural units
double L;

//  Initial Temperature in Natural Units
double Tinit;  //2;

//  Vectors!
//  Position
double* r = (double *) malloc(VECSIZE*sizeof(double));
//  Velocity
double* v= (double *) malloc(VECSIZE*sizeof(double));
//  Acceleration
double* a= (double *) malloc(VECSIZE*sizeof(double));


double *d_r, *d_a; 
double aux = N * 3 * sizeof(double);

char *atype = (char *)malloc(3 * sizeof(char));


double gaussdist() {
    static bool available = false;
    static double gset;
    double fac, rsq, v1, v2, returnValue;
    if (!available) {
        do {
            v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        available = true;
        
        returnValue =  v2*fac;
    } else {
        
        available = false;
        returnValue = gset;
        
    }
    return returnValue;
}

void initializeVelocities() {
    int i;
    double vCM[3] = {0, 0, 0};
    for (i=0; i < VECSIZE;i += 3) {
        v[i] = gaussdist();
        v[i+1] = gaussdist();
        v[i+2] = gaussdist();
        vCM[0] += v[i];
        vCM[1] += v[i+1];
        vCM[2] += v[i+2];
    }
    vCM[0] /= N;
    vCM[1] /= N;
    vCM[2] /= N;
   
    //  Now we want to scale the average velocity of the system
    //  by a factor which is consistent with our initial temperature, Tinit
    double vSqdSum, lambda;
    vSqdSum=0.;
    for (i=0; i < VECSIZE;i += 3) {
        v[i] -= vCM[0];
        v[i+1] -= vCM[1];
        v[i+2] -= vCM[2];
        vSqdSum += (v[i]*v[i] + v[i + 1]*v[i + 1] + v[i + 2]*v[i + 2]);
    }
    
    lambda = sqrt(3*(N-1)*Tinit/vSqdSum);
    
    for (i=0; i < VECSIZE;) {
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
    double velo_1 = 0., velo_2 = 0.;
    for(int i=0; i < VECSIZE; i+=2){
        velo_1 += v[i]*v[i]; 
        velo_2 += v[i+1]*v[i+1];
    }
    double velo = velo_1 + velo_2;
    KE = velo/2;
    mvs = velo/N;
}

__device__ double atomicAddDouble(double* address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;

    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed, __double_as_longlong(val +
                       __longlong_as_double(assumed)));
    } while (assumed != old);

    return __longlong_as_double(old);
}

 
__global__ void computeAccelerations_plus_potentialGPU(double *d_a, double *d_r, double *d_Pot) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    __shared__ double sharedRk[NUM_THREADS_PER_BLOCK * 3];

    sharedRk[threadIdx.x * 3] = d_r[i * 3];
    sharedRk[threadIdx.x * 3 + 1] = d_r[i * 3 + 1];
    sharedRk[threadIdx.x * 3 + 1] = d_r[i * 3 + 2];
    

    if (i < NCuda - 1) {
        double vPot_local = 0.0;
        double a0 = 0.0, a1 = 0.0, a2 = 0.0;

        for (int j = i + 1; j < NCuda; j++) {
            double M0 = sharedRk[threadIdx.x * 3] - d_r[j * 3],
                   M1 = sharedRk[threadIdx.x * 3 + 1] - d_r[j * 3 + 1],
                   M2 = sharedRk[threadIdx.x * 3 + 2] - d_r[j * 3 + 2];

            double rSqd = M0 * M0 + M1 * M1 + M2 * M2;

            double rSqd3 = rSqd * rSqd * rSqd;
            vPot_local += (1 - rSqd3) / (rSqd3 * rSqd3);
            double f = (48. - 24. * rSqd3) / (rSqd3 * rSqd3 * rSqd);

            double aux0 = M0 * f;
            double aux1 = M1 * f;
            double aux2 = M2 * f;

            a0 += aux0;
            a1 += aux1;
            a2 += aux2;

            atomicAddDouble(&d_a[j * 3], -aux0);
            atomicAddDouble(&d_a[j * 3 + 1], -aux1);
            atomicAddDouble(&d_a[j * 3 + 2], -aux2);

        }
        d_Pot[i] = vPot_local;

        atomicAddDouble(&d_a[i * 3], a0);
        atomicAddDouble(&d_a[i * 3 + 1], a1);
        atomicAddDouble(&d_a[i * 3 + 2], a2);
    }
}

void computeAccelerations_plus_potential() {

    double Pot = 0.0;
    double v_Pot[N];
    double* d_Pot;

    cudaMalloc((void**)&d_r, aux);
    cudaMalloc((void**)&d_a, aux);
    cudaMalloc((void**)&d_Pot, N * sizeof(double) - 1);
    checkCUDAError("Mem Allocation");

    cudaMemcpy(d_a, a, aux, cudaMemcpyHostToDevice);
    cudaMemcpy(d_r, r, aux, cudaMemcpyHostToDevice);
    checkCUDAError("Memcpy Host -> Device");

    int bpg = (N + NUM_THREADS_PER_BLOCK - 1) / NUM_THREADS_PER_BLOCK; 
    computeAccelerations_plus_potentialGPU<<<bpg, NUM_THREADS_PER_BLOCK>>>(d_a, d_r, d_Pot);
    cudaDeviceSynchronize();
    checkCUDAError("Error in PotentialGPU");

    cudaMemcpy(a, d_a, aux, cudaMemcpyDeviceToHost);

    cudaMemcpy(v_Pot, d_Pot, N * sizeof(double) - 1, cudaMemcpyDeviceToHost);
    checkCUDAError("Memcpy Device -> Host");

    for (int i = 0; i < N; i++) {
        Pot += v_Pot[i];
    }

    cudaFree(d_r);
    cudaFree(d_a);
    cudaFree(d_Pot);
    checkCUDAError("Free Mem");

    PE = Pot * 8;
}


// returns sum of dv/dt*m/A (aka Pressure) from elastic collisions with walls
double VelocityVerlet(double dt) {
    double psum = 0.0, half_dt = 0.5*dt;
    //  Compute accelerations from forces at current position
    // this call was removed (commented) for prledagogical reasons
    //  Update positions and velocity with current velocity and acceleration
    for (int i=0; i < VECSIZE; i+=3) {
        v[i] += a[i] * half_dt;
        v[i+1] += a[i+1] * half_dt;
        v[i+2] += a[i+2] * half_dt;
        r[i] += v[i] * dt;     
        r[i+1] += v[i+1] * dt; 
        r[i+2] += v[i+2] * dt; 
        a[i] = 0.0;
        a[i+1] = 0.0;
        a[i+2] = 0.0;
    }

    computeAccelerations_plus_potential();

    //  Update velocity with updated acceleration
    for (int i=0; i < VECSIZE; i++){
        v[i] += a[i] * half_dt;
        // Elastic walls
        if (r[i]<0. || r[i]>=L) {
            v[i] *=-1.; //- elastic walls
            psum += fabs(v[i]);  // contribution to pressure from "left" walls
        }
    }
    return psum/(3*L*L*dt);
}

int main(){
    int i;
    double Vol, Temp, Press, rho;
    double VolFac, TempFac, PressFac, timefac;
    char prefix[1000], tfn[1000], ofn[1000], afn[1000];
    FILE *tfp, *ofp, *afp;
    
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("                  WELCOME TO WILLY P CHEM MD!\n");
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n  ENTER A TITLE FOR YOUR CALCULATION!\n");
    scanf("%s",prefix);
    strcpy(tfn,prefix);
    strcat(tfn,"_traj.xyz");
    strcpy(ofn,prefix);
    strcat(ofn,"_output.txt");
    strcpy(afn,prefix);
    strcat(afn,"_average.txt");
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("                  TITLE ENTERED AS '%s'\n",prefix);
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("  WHICH NOBLE GAS WOULD YOU LIKE TO SIMULATE? (DEFAULT IS ARGON)\n");
    printf("\n  FOR HELIUM,  TYPE 'He' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR NEON,    TYPE 'Ne' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR ARGON,   TYPE 'Ar' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR KRYPTON, TYPE 'Kr' THEN PRESS 'return' TO CONTINUE\n");
    printf("  FOR XENON,   TYPE 'Xe' THEN PRESS 'return' TO CONTINUE\n");
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    scanf("%s",atype);
    
    if (strcmp(atype,"He")==0) {
        VolFac = 1.8399744000000005e-29;
        PressFac = 8152287.336171632;
        TempFac = 10.864459551225972;
        timefac = 1.7572698825166272e-12;
    }else if (strcmp(atype,"Ne")==0) {  
        VolFac = 2.0570823999999997e-29;
        PressFac = 27223022.27659913;
        TempFac = 40.560648991243625;
        timefac = 2.1192341945685407e-12;  
    }else if (strcmp(atype,"Ar")==0) {
        VolFac = 3.7949992920124995e-29;
        PressFac = 51695201.06691862;
        TempFac = 142.0950000000000;
        timefac = 2.09618e-12;  
    }else if (strcmp(atype,"Kr")==0) {
        VolFac = 4.5882712000000004e-29;
        PressFac = 59935428.40275003;
        TempFac = 199.1817584391428;
        timefac = 8.051563913585078e-13;
    }else if (strcmp(atype,"Xe")==0) {
        VolFac = 5.4872e-29;
        PressFac = 70527773.72794868;
        TempFac = 280.30305642163006;
        timefac = 9.018957925790732e-13;   
    }else {
        VolFac = 3.7949992920124995e-29;
        PressFac = 51695201.06691862;
        TempFac = 142.0950000000000;
        timefac = 2.09618e-12;
        strcpy(atype,"Ar");   
    }
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n                     YOU ARE SIMULATING %s GAS! \n",atype);
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n  YOU WILL NOW ENTER A FEW SIMULATION PARAMETERS\n");
    printf("  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
    printf("\n\n  ENTER THE INTIAL TEMPERATURE OF YOUR GAS IN KELVIN\n");
    scanf("%lf",&Tinit);

    if (Tinit<0.) {
        printf("\n  !!!!! ABSOLUTE TEMPERATURE MUST BE A POSITIVE NUMBER!  PLEASE TRY AGAIN WITH A POSITIVE TEMPERATURE!!!\n");
        exit(0);
    }

    if (N>=VECSIZE) {
        printf("\n\n\n  MAXIMUM NUMBER OF PARTICLES IS %i\n\n  PLEASE ADJUST YOUR INPUT FILE ACCORDINGLY \n\n", VECSIZE);
        exit(0);
    }

    // Convert initial temperature from kelvin to natural units
    Tinit /= TempFac;
    scanf("%lf",&rho);

    // Copy N to the device variable d_NCuda
    cudaMemcpyToSymbol(NCuda, &N, sizeof(int));

    Vol = N/(rho*NA);
    Vol /= VolFac;

    if (Vol<N) {
        printf("\n\n\n  YOUR DENSITY IS VERY HIGH!\n\n");
        printf("  THE NUMBER OF PARTICLES IS %i AND THE AVAILABLE VOLUME IS %f NATURAL UNITS\n",N,Vol);
        printf("  SIMULATIONS WITH DENSITY GREATER THAN 1 PARTCICLE/(1 Natural Unit of Volume) MAY DIVERGE\n");
        printf("  PLEASE ADJUST YOUR INPUT FILE ACCORDINGLY AND RETRY\n\n");
        exit(0);
    }

    // Length of the box in natural units:
    L = pow(Vol,(1./3));
    
    //  Files that we can write different quantities to
    tfp = fopen(tfn,"w");    //  The MD trajectory, coordinates of every particle at each timestep
    ofp = fopen(ofn,"w");    //  Output of other quantities (T, P, gc, etc) at every timestep
    afp = fopen(afn,"w");    //  Average T, P, gc, etc from the simulation
    
    int NumTime = 200;
    double dt = 0.5e-14/timefac;

    if (strcmp(atype,"He")==0) {
        dt = 0.2e-14/timefac;
        NumTime=50000;
    }

    initialize();
    
    computeAccelerations_plus_potential();

    fprintf(tfp,"%i\n",N);

    double Pavg = 0.0;
    double Tavg = 0.0;
    fprintf(ofp,"  time (s)              T(t) (K)              P(t) (Pa)           Kinetic En. (n.u.)     Potential En. (n.u.) Total En. (n.u.)\n");
    for (i=0; i<NumTime+1; i++) {
        //  This just prints updates on progress of the calculation for the users convenience
        
        Press = VelocityVerlet(dt)* PressFac;

        MeanSquaredVelocity_and_Kinetic();

        Temp = mvs/3 * TempFac;

        Tavg += Temp;
        Pavg += Press;
        
        fprintf(ofp,"  %8.4e  %20.8f  %20.8f %20.8f  %20.8f  %20.8f \n",i*dt*timefac,Temp,Press,KE, PE, KE+PE);
    }
    
    // Because we have calculated the instantaneous temperature and pressure,
    // we can take the average over the whole simulation here
    Pavg /= NumTime;
    Tavg /= NumTime;
    double Z = Pavg*(Vol*VolFac)/(N*kBSI*Tavg);
    double gc = NA*Pavg*(Vol*VolFac)/(N*Tavg);
    fprintf(afp,"  Total Time (s)      T (K)               P (Pa)      PV/nT (J/(mol K))         Z           V (m^3)              N\n");
    fprintf(afp," --------------   -----------        ---------------   --------------   ---------------   ------------   -----------\n");
    fprintf(afp,"  %8.4e  %15.5f       %15.5f     %10.5f       %10.5f        %10.5e         %i\n",i*dt*timefac,Tavg,Pavg,gc,Z,Vol*VolFac,N);
    
    printf("\n  TO ANIMATE YOUR SIMULATION, OPEN THE FILE \n  '%s' WITH VMD AFTER THE SIMULATION COMPLETES\n",tfn);
    printf("\n  TO ANALYZE INSTANTANEOUS DATA ABOUT YOUR MOLECULE, OPEN THE FILE \n  '%s' WITH YOUR FAVORITE TEXT EDITOR OR IMPORT THE DATA INTO EXCEL\n",ofn);
    printf("\n  THE FOLLOWING THERMODYNAMIC AVERAGES WILL BE COMPUTED AND WRITTEN TO THE FILE  \n  '%s':\n",afn);
    printf("\n  AVERAGE TEMPERATURE (K):                 %15.5f\n",Tavg);
    printf("\n  AVERAGE PRESSURE  (Pa):                  %15.5f\n",Pavg);
    printf("\n  PV/nT (J * mol^-1 K^-1):                 %15.5f\n",gc);
    printf("\n  PERCENT ERROR of pV/nT AND GAS CONSTANT: %15.5f\n",100*fabs(gc-8.3144598)/8.3144598);
    printf("\n  THE COMPRESSIBILITY (unitless):          %15.5f \n",Z);
    printf("\n  TOTAL VOLUME (m^3):                      %10.5e \n",Vol*VolFac);
    printf("\n  NUMBER OF PARTICLES (unitless):          %i \n", N); 
    
    free(r);
    free(v);
    free(a);
    free(atype);
    fclose(tfp);
    fclose(ofp);
    fclose(afp);
    
    return 0;
}