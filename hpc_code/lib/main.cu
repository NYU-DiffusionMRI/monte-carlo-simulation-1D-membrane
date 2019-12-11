//
//  main.cpp
//  diffusion_membrane_exchange
//
//  Update Journal:
//  -- 11/14/2019: re-write cylinder code to membrane code
//  -- 11/17/2019: membrane thickness = 0
//
//  Created by Hong-Hsi Lee in November, 2019.
//


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <iomanip>
#include <time.h>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <complex>
#include <string>
//#include "diffusion_lib.h"

#include <cuda.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/copy.h>

using namespace std;
    
#define Pi 3.14159265
#define timepoints 10000

// ********** diffusion kernel **********
__device__ void pixPosition ( const double &x_in, const unsigned int &NPix, int &xPix ) {
    double x=x_in;
    if ( x<0 ) { x+=1; }
    if ( x>1 ) { x-=1; }
    xPix=floor(x*NPix);
}

__device__ void translateXc ( const double &x, double &xc ) {
    // Translate circle center xc to make it as close to the position x as possible
    int ii=2*(xc<0.5)-1;
    double ti=0, d2=fabs(x-xc), d2Tmp=fabs(x-xc-ii);
    if (d2Tmp<d2) { ti=ii; }
    xc+=ti;
}

__device__ void crossmembrane( const double &xi, const double &xf, const double &xc_in, const bool &translateFlag, bool &instruction) {
    double xc=xc_in;
    if (translateFlag) { translateXc(xf,xc); }
    instruction = (xc-xi)*(xc-xf)<=0;
}

__device__ void elastic(const double &x, const double &v, const double &dx, const double &xc_in, const bool &translateFlag, double &xt) {
    double xc=xc_in;
    if (translateFlag) { translateXc(x,xc); }
    if (fabs(x-xc)<dx) {
        xt=xc-v*(dx-fabs(x-xc));
    }
}

// ********** cuda kernel **********
__device__ double atomAdd(double* address, double val)
{
    unsigned long long int* address_as_ull =
    (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                                             __longlong_as_double(assumed)));
        
        // Note: uses integer comparison to avoid hang in case of NaN (since NaN != NaN)
    } while (assumed != old);
    
    return __longlong_as_double(old);
}

__global__ void setup_kernel(curandStatePhilox4_32_10_t *state, unsigned long seed){
    int idx = threadIdx.x+blockDim.x*blockIdx.x;
    curand_init(seed, idx, 0, &state[idx]);
}

__global__ void propagate(curandStatePhilox4_32_10_t *state, double *dx2, double *dx4, double *sig, const int TN, const int NPar, const double res, const double dt, const double stepIN, const double probI, const unsigned int NPix, const int initFlag, const double *xCir, const bool *translateFlag, const int Nbval, const double *bval, const unsigned int *APix){
    int idx = threadIdx.x + blockDim.x * blockIdx.x;
    int stride = blockDim.x * gridDim.x;
    curandStatePhilox4_32_10_t localstate=state[idx];
    
    int Tstep=TN/timepoints;
    
    for (int k=idx; k<NPar; k+=stride){
        double xPar=0, xCirTmp=0;
        
        unsigned int a=0, aTmp=0;                           // Element of APix matrix
        unsigned int acell[2]={0}; bool instruction[2]={0}; // Cell label
        
        double xi=0, xt=0, xTmp=0;                          // Particle position
        int xtG=0, xTmpG=0;                                 // Position on grid after diffusion
        double vrand=0;                                     // Random number
        int tidx=0;
        
        double vp=0;                                        // Nomalized diffusion velocity
        int acell_hit=0;                                    // Label of the cell encountered by the walker
        
        double xjmp=0;
        double dx=0, TD=0, qx=0;
        
        //********** Initialize Walker Positions *********
        while (1){
            xPar=curand_uniform_double(&localstate);
            
            if ( initFlag==1 ) { // 1. Initial positon: ICS
                xi=xPar;
                break;
            } else if ( initFlag==2 ) { // 2. Initial position: center
                xi=0.5;
                break;
            }
        }
        
        // ********** Simulate diffusion **********
        xt=xi;
        for (int i=0; i<TN; i++){
            // The cells close to the walker in the previous step
            pixPosition(xt,NPix,xtG);                             // Position on grid
            a=APix[ xtG ];
            
            // ********** One step **********
            // Case 1 Diffusion In ICS
            acell[0]=a;
            // Primitive position after diffusion
            vrand=curand_uniform_double(&localstate);
            vp=2*(vrand>0.5)-1;
            xTmp=xt+stepIN*vp;
            
            pixPosition(xTmp,NPix,xTmpG);
            aTmp=APix[ xTmpG ];
            acell[1]=aTmp;
            
            // Check if the segment(xt,xTmp) overlaps with any membrane
            for (int j=0; j<2; j++) {
                if (acell[j]==0) {
                    instruction[j]=false;
                } else {
                    crossmembrane(xt,xTmp,xCir[acell[j]-1],translateFlag[acell[j]-1],instruction[j]);
                }
            }
            if ( (instruction[0]==false) & (instruction[1]==false) ) {
                // Case 1.1 Walker diffuses in ICS and does not encounter any cell membrane.
                xt=xTmp;
            } else {
                // Case 1.2 Walker diffuses in ICS and encounters the cell membrane.
                
                // Determine the cell to collide with.
                acell_hit=0;
                for (int j=0; j<2; j++) {
                    if ( instruction[j] ) { acell_hit=acell[j]; }
                }
                
                if ( acell_hit==0 ){
                    printf("error: Walker in ICS does not encounter the cell membrane.\n");
                }
                
                xCirTmp=xCir[acell_hit-1];
                
                vrand=curand_uniform_double(&localstate);
                if (vrand<probI) {
                    // Case 1.2.1 Permeation through membrane
                    xt=xTmp;
                }
                else {
                    // Case 1.2.2 Elastic collision in ECS
                    elastic(xt, vp, stepIN, xCirTmp, translateFlag[acell_hit-1], xTmp);
                    xt=xTmp;
                }
            }
            
            // Periodic boundary condition
            if (xt>1) {
                xt-=1;
                xjmp+=1;
            }
            else if (xt<0) {
                xt+=1;
                xjmp-=1;
            }
            
            // ********** End one step **********
            
            if ( (i%Tstep) ==0 ) { // Save moment tensor for dx^2 and dx^4, and signal for the b-table
                tidx=i/Tstep;
                
                dx=(xt+xjmp-xi)*res;
                atomAdd(&dx2[tidx],dx*dx);
                atomAdd(&dx4[tidx],dx*dx*dx*dx);
                
                TD=i*dt;
                for (int j=0; j<Nbval; j++) {
                    qx = sqrt(bval[j]/TD) * dx;
                    atomAdd(&sig[tidx*Nbval+j],cos(qx));
                }
            }
            
        }
    }
    state[idx]=localstate;
}

    
//********** Define tissue parameters **********

int main(int argc, char *argv[]) {
    
    clock_t begin=clock();
    clock_t end=clock();
    
    // Define index number
    int i=0, j=0;
    
    //********** Load mictostructure **********
    
    double dt=0;                // Time step in ms
    int TN=0;                   // Number of time steps
    int NPar=0;                 // Number of particles
    int Nbval=0;                // Number of b-values
    
    double Din=0;               // Diffusion coefficient inside the axon in µm^2/ms
    double kappa=0;             // Permeability of a lipid bi-layer in µm/ms
    int initFlag=1;             // Initial position: 1=ICS, 2=center
    int thread_per_block=0;     // Number of threads per block
    
    unsigned int NPix=0, NAx=0;
    double res=0;
    
    // simulation parameter
    ifstream myfile0 ("simParamInput.txt", ios::in);
    myfile0>>dt; myfile0>>TN; myfile0>>NPar;
    myfile0>>Din; myfile0>>kappa;
    myfile0>>initFlag;
    myfile0>>thread_per_block;
    myfile0.close();
    
    double stepIN=sqrt(2.0*dt*Din);     // Step size in ICS in µm
    
    // resolution
    ifstream myfile1 ("phantom_res.txt", ios::in);
    myfile1>>res;
    myfile1.close();
    
    // Pixel # along each side
    ifstream myfile2 ("phantom_NPix.txt", ios::in);
    myfile2>>NPix;
    myfile2.close();
    
    // Pixelized matrix A indicating axon labels
    thrust::host_vector<unsigned int> APix(NPix);
    ifstream myfile3 ("phantom_APix.txt", ios::in);
    for (i=0; i<NPix; i++){
        myfile3>>APix[i];
    }
    myfile3.close();
    
    // Number of axons
    ifstream myfile4 ("phantom_NMem.txt", ios::in);
    myfile4>>NAx;
    myfile4.close();
    
    // Circle center of x-coordinate
    thrust::host_vector<double> xCir(NAx);
    ifstream myfile5 ("phantom_xMem.txt", ios::in);
    for (i=0; i<NAx; i++){
        myfile5>>xCir[i];
    }
    myfile5.close();
    
    // Number of b-values
    ifstream myfile6 ("Nbval.txt", ios::in);
    myfile6>>Nbval;
    myfile6.close();
    
    // b-value
    thrust::host_vector<double> bval(Nbval);
    ifstream myfile7 ("bval.txt", ios::in);
    for (i=0; i<Nbval; i++){
        myfile7>>bval[i];
    }
    myfile7.close();
    
    //********** Initialize Particle Positions in IAS *********
    const double probI=stepIN*kappa/Din;          // Probability constant from ICS to myelin
    stepIN/=res;
    cout<<"probI="<<probI<<endl;
    cout<<"NPix="<<NPix<<endl;
    cout<<"NAx="<<NAx<<endl;
    // Create translate flag to speed up the code
    thrust::host_vector<bool> translateFlag(NAx);
    for (i=0; i<NAx; i++) {
        if ( (xCir[i]+2*stepIN>=1) | (xCir[i]-2*stepIN<=0) ) {
            translateFlag[i]=true;
        } else {
            translateFlag[i]=false;
        }
    }
    
    // ********** Simulate diffusion **********
    
    // Initialize seed
    unsigned long seed=0;
    FILE *urandom;
    urandom = fopen("/dev/random", "r");
    fread(&seed, sizeof (seed), 1, urandom);
    fclose(urandom);
    
    // Initialize state of RNG
    int blockSize = thread_per_block;
    int numBlocks = (NPar + blockSize - 1) / blockSize;
    cout<<numBlocks<<endl<<blockSize<<endl;
    
    thrust::device_vector<curandStatePhilox4_32_10_t> devState(numBlocks*blockSize);
    setup_kernel<<<numBlocks, blockSize>>>(devState.data().get(),seed);
    
    // Initialize output
    thrust::host_vector<double> dx2(timepoints);
    thrust::host_vector<double> dx4(timepoints);
    thrust::host_vector<double> sig(timepoints*Nbval);
    for (i=0;i<timepoints;i++){ dx2[i]=0; }
    for (i=0;i<timepoints;i++){ dx4[i]=0; }
    for (i=0;i<timepoints*Nbval;i++){ sig[i]=0; }
    
    // Move data from host to device
    thrust::device_vector<double> d_dx2=dx2;
    thrust::device_vector<double> d_dx4=dx4;
    thrust::device_vector<double> d_sig=sig;
    thrust::device_vector<double> d_xCir=xCir;
    thrust::device_vector<bool> d_translateFlag=translateFlag;
    thrust::device_vector<double> d_bval=bval;
    thrust::device_vector<unsigned int> d_APix=APix;
    
    // Parallel computation
    begin=clock();
    propagate<<<numBlocks, blockSize>>>(devState.data().get(), d_dx2.data().get(), d_dx4.data().get(), d_sig.data().get(), TN, NPar, res, dt, stepIN, probI, NPix, initFlag, d_xCir.data().get(), d_translateFlag.data().get(), Nbval, d_bval.data().get(), d_APix.data().get());
    cudaDeviceSynchronize();
    end=clock();
    cout << "Done! Elpased time "<<double((end-begin)/CLOCKS_PER_SEC) << " s"<< endl;
    
    thrust::copy(d_dx2.begin(), d_dx2.end(), dx2.begin());
    thrust::copy(d_dx4.begin(), d_dx4.end(), dx4.begin());
    thrust::copy(d_sig.begin(), d_sig.end(), sig.begin());
    
    ofstream fdx2out("dx2_diffusion.txt");
    ofstream fdx4out("dx4_diffusion.txt");
    ofstream fsigout("sig_diffusion.txt");
    fdx2out.precision(15);
    fdx4out.precision(15);
    fsigout.precision(15);
    
    for (i=0; i<timepoints; i++) {
        fdx2out<<dx2[i]<<endl;
        fdx4out<<dx4[i]<<endl;
        for (j=0; j<Nbval; j++) {
            if ( j==(Nbval-1) ) {
                fsigout<<sig[i*Nbval+j]<<endl;
            } else {
                fsigout<<sig[i*Nbval+j]<<"\t";
            }
        }
    }
    fdx2out.close();
    fdx4out.close();
    fsigout.close();
    
    ofstream paraout ("sim_para.txt");
    paraout<<dt<<endl<<TN<<endl<<NPar<<endl;
    paraout<<Din<<endl;
    paraout<<kappa<<endl<<initFlag<<endl;
    paraout.close();
    
    ofstream TDout("diff_time.txt");
    for (i=0; i<timepoints; i++){
        TDout<<(i*(TN/timepoints)+1)*dt<<endl;
    }
    TDout.close();
}

