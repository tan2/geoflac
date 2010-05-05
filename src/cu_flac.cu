// -*- C++ -*-

#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>

#include "cu_flac.h"

// size of thread block
// requirement: nthz is multiple of 8, nthx*nthz is multiple of 32
static const int nthx = 16;
static const int nthz = 16;

// global variable holding constant model parameters
flac_param_t param;


/*
 * Error reporting for CUDA calls
 */
__host__ static
void checkCUDAError(const char *msg)
{
    cudaError_t err;

    // uncomment the following line to debug cuda calls ...
    //cudaThreadSynchronize();

    err = cudaGetLastError();
    if(cudaSuccess != err) {
        fprintf(stderr, "CUDA error: %s: %s.\n", msg,
                cudaGetErrorString(err));
        exit(99);
    }
    return;
}


__global__ static
void cu_fl_node(double *force_d, double *balance_d, double *vel_d,
                const double *cord_d, const double *stress0_d, const double *temp_d,
                const double *rmass_d, const double *amass_d,
                const double *bc_d, const int *ncod_d,
                const double dt, const double drat, const double time,
                int nx, int nz)
{
    // i and j are indices to fortran arrays
    const int i = blockIdx.x * blockDim.x + threadIdx.x + 1;
    const int j = blockIdx.y * blockDim.y + threadIdx.y + 1;

    double fx, fy;
    double fcx, fcy, blx, bly;

    if(j > nz || i > nx) return;

    if(PAR.ynstressbc) {
        fcx = force(j,i,1);
        fcy = force(j,i,2);
        blx = balance(j,i,1);
        bly = balance(j,i,2);
    } else {
        fcx = fcy = blx = bly = 0.0;
    }

    // REGULAR PART - forces from stresses

    // Element (j-1,i-1). Triangles B,C,D
    if (j>1 && i>1 ) {

        // triangle B
        // side 2-3
        fx = stress0(j-1,i-1,1,2) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) -
            stress0(j-1,i-1,3,2) * (cord(j  ,i  ,1)-cord(j  ,i-1,1));
        fy = stress0(j-1,i-1,3,2) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) -
            stress0(j-1,i-1,2,2) * (cord(j  ,i  ,1)-cord(j  ,i-1,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
        // side 3-1
        fx = stress0(j-1,i-1,1,2) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) -
            stress0(j-1,i-1,3,2) * (cord(j-1,i  ,1)-cord(j  ,i  ,1));
        fy = stress0(j-1,i-1,3,2) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) -
            stress0(j-1,i-1,2,2) * (cord(j-1,i  ,1)-cord(j  ,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);

        // triangle C
        // side 2-3
        fx = stress0(j-1,i-1,1,3) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) -
            stress0(j-1,i-1,3,3) * (cord(j  ,i  ,1)-cord(j  ,i-1,1));
        fy = stress0(j-1,i-1,3,3) * (cord(j  ,i  ,2)-cord(j  ,i-1,2)) -
            stress0(j-1,i-1,2,3) * (cord(j  ,i  ,1)-cord(j  ,i-1,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
        // side 3-1
        fx = stress0(j-1,i-1,1,3) * (cord(j-1,i-1,2)-cord(j  ,i  ,2)) -
         stress0(j-1,i-1,3,3) * (cord(j-1,i-1,1)-cord(j  ,i  ,1));
        fy = stress0(j-1,i-1,3,3) * (cord(j-1,i-1,2)-cord(j  ,i  ,2)) -
            stress0(j-1,i-1,2,3) * (cord(j-1,i-1,1)-cord(j  ,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);

        // triangle D
        // side 1-2
        fx = stress0(j-1,i-1,1,4) * (cord(j  ,i  ,2)-cord(j-1,i-1,2)) -
            stress0(j-1,i-1,3,4) * (cord(j  ,i  ,1)-cord(j-1,i-1,1));
        fy = stress0(j-1,i-1,3,4) * (cord(j  ,i  ,2)-cord(j-1,i-1,2)) -
            stress0(j-1,i-1,2,4) * (cord(j  ,i  ,1)-cord(j-1,i-1,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
        // side 2-3
        fx = stress0(j-1,i-1,1,4) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) -
            stress0(j-1,i-1,3,4) * (cord(j-1,i  ,1)-cord(j  ,i  ,1));
        fy = stress0(j-1,i-1,3,4) * (cord(j-1,i  ,2)-cord(j  ,i  ,2)) -
            stress0(j-1,i-1,2,4) * (cord(j-1,i  ,1)-cord(j  ,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
    }

    // Element (j-1,i). Triangles A,B,C.
    if (j>1 && i<nx) {

        // triangle A
        // side 1-2
        fx = stress0(j-1,i  ,1,1) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) -
            stress0(j-1,i  ,3,1) * (cord(j  ,i  ,1)-cord(j-1,i  ,1));
        fy = stress0(j-1,i  ,3,1) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) -
            stress0(j-1,i  ,2,1) * (cord(j  ,i  ,1)-cord(j-1,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
        // side 2-3
        fx = stress0(j-1,i  ,1,1) * (cord(j-1,i+1,2)-cord(j  ,i  ,2)) -
            stress0(j-1,i  ,3,1) * (cord(j-1,i+1,1)-cord(j  ,i  ,1));
        fy = stress0(j-1,i  ,3,1) * (cord(j-1,i+1,2)-cord(j  ,i  ,2)) -
            stress0(j-1,i  ,2,1) * (cord(j-1,i+1,1)-cord(j  ,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);

        // triangle B
        // side 1-2
        fx = stress0(j-1,i  ,1,2) * (cord(j  ,i  ,2)-cord(j-1,i+1,2)) -
            stress0(j-1,i  ,3,2) * (cord(j  ,i  ,1)-cord(j-1,i+1,1));
        fy = stress0(j-1,i  ,3,2) * (cord(j  ,i  ,2)-cord(j-1,i+1,2)) -
            stress0(j-1,i  ,2,2) * (cord(j  ,i  ,1)-cord(j-1,i+1,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
        // side 2-3
        fx = stress0(j-1,i  ,1,2) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) -
            stress0(j-1,i  ,3,2) * (cord(j  ,i+1,1)-cord(j  ,i  ,1));
        fy = stress0(j-1,i  ,3,2) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) -
            stress0(j-1,i  ,2,2) * (cord(j  ,i+1,1)-cord(j  ,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);

        // triangle C
        // side 1-2
        fx = stress0(j-1,i  ,1,3) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) -
            stress0(j-1,i  ,3,3) * (cord(j  ,i  ,1)-cord(j-1,i  ,1));
        fy = stress0(j-1,i  ,3,3) * (cord(j  ,i  ,2)-cord(j-1,i  ,2)) -
            stress0(j-1,i  ,2,3) * (cord(j  ,i  ,1)-cord(j-1,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
        // side 2-3
        fx = stress0(j-1,i  ,1,3) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) -
            stress0(j-1,i  ,3,3) * (cord(j  ,i+1,1)-cord(j  ,i  ,1));
        fy = stress0(j-1,i  ,3,3) * (cord(j  ,i+1,2)-cord(j  ,i  ,2)) -
            stress0(j-1,i  ,2,3) * (cord(j  ,i+1,1)-cord(j  ,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
    }

    // Element (j,i-1). Triangles A,B,D
    if (j<nz && i>1) {

        // triangle A
        // side 2-3
        fx = stress0(j  ,i-1,1,1) * (cord(j  ,i  ,2)-cord(j+1,i-1,2)) -
            stress0(j  ,i-1,3,1) * (cord(j  ,i  ,1)-cord(j+1,i-1,1));
        fy = stress0(j  ,i-1,3,1) * (cord(j  ,i  ,2)-cord(j+1,i-1,2)) -
            stress0(j  ,i-1,2,1) * (cord(j  ,i  ,1)-cord(j+1,i-1,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
        // side 3-1
        fx = stress0(j  ,i-1,1,1) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) -
            stress0(j  ,i-1,3,1) * (cord(j  ,i-1,1)-cord(j  ,i  ,1));
        fy = stress0(j  ,i-1,3,1) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) -
            stress0(j  ,i-1,2,1) * (cord(j  ,i-1,1)-cord(j  ,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);

        // triangle B
        // side 1-2
        fx = stress0(j  ,i-1,1,2) * (cord(j+1,i-1,2)-cord(j  ,i  ,2)) -
            stress0(j  ,i-1,3,2) * (cord(j+1,i-1,1)-cord(j  ,i  ,1));
        fy = stress0(j  ,i-1,3,2) * (cord(j+1,i-1,2)-cord(j  ,i  ,2)) -
            stress0(j  ,i-1,2,2) * (cord(j+1,i-1,1)-cord(j  ,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
        // side 3-1
        fx = stress0(j  ,i-1,1,2) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) -
            stress0(j  ,i-1,3,2) * (cord(j  ,i  ,1)-cord(j+1,i  ,1));
        fy = stress0(j  ,i-1,3,2) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) -
            stress0(j  ,i-1,2,2) * (cord(j  ,i  ,1)-cord(j+1,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);

        // triangle D
        // side 2-3
        fx = stress0(j  ,i-1,1,4) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) -
            stress0(j  ,i-1,3,4) * (cord(j  ,i  ,1)-cord(j+1,i  ,1));
        fy = stress0(j  ,i-1,3,4) * (cord(j  ,i  ,2)-cord(j+1,i  ,2)) -
            stress0(j  ,i-1,2,4) * (cord(j  ,i  ,1)-cord(j+1,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
        // side 3-1
        fx = stress0(j  ,i-1,1,4) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) -
            stress0(j  ,i-1,3,4) * (cord(j  ,i-1,1)-cord(j  ,i  ,1));
        fy = stress0(j  ,i-1,3,4) * (cord(j  ,i-1,2)-cord(j  ,i  ,2)) -
            stress0(j  ,i-1,2,4) * (cord(j  ,i-1,1)-cord(j  ,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
    }

    // Element (j,i). Triangles A,C,D
    if (j<nz && i<nx ) {

        // triangle A
        // side 1-2
        fx = stress0(j  ,i  ,1,1) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) -
            stress0(j  ,i  ,3,1) * (cord(j+1,i  ,1)-cord(j  ,i  ,1));
        fy = stress0(j  ,i  ,3,1) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) -
            stress0(j  ,i  ,2,1) * (cord(j+1,i  ,1)-cord(j  ,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
        // side 3-1
        fx = stress0(j  ,i  ,1,1) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) -
            stress0(j  ,i  ,3,1) * (cord(j  ,i  ,1)-cord(j  ,i+1,1));
        fy = stress0(j  ,i  ,3,1) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) -
            stress0(j  ,i  ,2,1) * (cord(j  ,i  ,1)-cord(j  ,i+1,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);

        // triangle C
        // side 1-2
        fx = stress0(j  ,i  ,1,3) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) -
            stress0(j  ,i  ,3,3) * (cord(j+1,i  ,1)-cord(j  ,i  ,1));
        fy = stress0(j  ,i  ,3,3) * (cord(j+1,i  ,2)-cord(j  ,i  ,2)) -
            stress0(j  ,i  ,2,3) * (cord(j+1,i  ,1)-cord(j  ,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
        // side 3-1
        fx = stress0(j  ,i  ,1,3) * (cord(j  ,i  ,2)-cord(j+1,i+1,2)) -
            stress0(j  ,i  ,3,3) * (cord(j  ,i  ,1)-cord(j+1,i+1,1));
        fy = stress0(j  ,i  ,3,3) * (cord(j  ,i  ,2)-cord(j+1,i+1,2)) -
            stress0(j  ,i  ,2,3) * (cord(j  ,i  ,1)-cord(j+1,i+1,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);

        // triangle D
        // side 1-2
        fx = stress0(j  ,i  ,1,4) * (cord(j+1,i+1,2)-cord(j  ,i  ,2)) -
            stress0(j  ,i  ,3,4) * (cord(j+1,i+1,1)-cord(j  ,i  ,1));
        fy = stress0(j  ,i  ,3,4) * (cord(j+1,i+1,2)-cord(j  ,i  ,2)) -
            stress0(j  ,i  ,2,4) * (cord(j+1,i+1,1)-cord(j  ,i  ,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
        // side 3-1
        fx = stress0(j  ,i  ,1,4) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) -
            stress0(j  ,i  ,3,4) * (cord(j  ,i  ,1)-cord(j  ,i+1,1));
        fy = stress0(j  ,i  ,3,4) * (cord(j  ,i  ,2)-cord(j  ,i+1,2)) -
            stress0(j  ,i  ,2,4) * (cord(j  ,i  ,1)-cord(j  ,i+1,1));
        fcx = fcx - 0.25*fx;
        fcy = fcy - 0.25*fy;
        blx = blx + 0.25*fabs(fx);
        bly = bly + 0.25*fabs(fy);
    }

    // GRAVITY FORCE
    fcy = fcy - rmass(j,i)*PAR.g;
    bly = bly + fabs(rmass(j,i)*PAR.g);


    if(PAR.nyhydro>0) {
        const int lneighbor = (i > 1) ? (i-1) : i;
        const int rneighbor = (i < nx) ? (i+1) : i;
        double dlx_l, dly_l, dlx_r, dly_r;
        double press_norm_l, press_norm_r;
        double p_est, rosubg;
        if((j == 1)) {
            // BOUNDARY CONDITIONS

            // pressure from water sea on top
            double rho_water = -10300.;
            double water_depth = 0.5*(cord(j,rneighbor,2)+cord(j,i,2));
            if (water_depth<0.) { // No water (above sea level)
                press_norm_l = rho_water*((cord(j,lneighbor,2)+cord(j,i,2))/2.);
                dlx_l = cord(j,i  ,1)-cord(j,lneighbor,1);
                dly_l = cord(j,i  ,2)-cord(j,lneighbor,2);

                press_norm_r = rho_water*((cord(j,rneighbor,2)+cord(j,i,2))/2.);
                dlx_r = cord(j,rneighbor,1)-cord(j,i,1);
                dly_r = cord(j,rneighbor,2)-cord(j,i,2);

                fcx = fcx - 0.5*press_norm_l*dly_l-0.5*press_norm_r*dly_r;
                fcy = fcy + 0.5*press_norm_l*dlx_l+0.5*press_norm_r*dlx_r;

                blx = 1.e+17;
            }
        }

        // bottom support - Archimed force (normal to the surface, shear component = 0)
        //write(*,*) force(nz,i,1),force(nz,i,2)
        if(j == nz) {
            p_est = PAR.pisos + 0.5*(den(PAR.iphsub)+PAR.drosub)*PAR.g*(cord(j,i,2)-PAR.rzbo);
            rosubg = PAR.g * (den(PAR.iphsub)+PAR.drosub) * (1-alfa(PAR.iphsub)*temp(j,i)+beta(PAR.iphsub)*p_est);

            press_norm_l = PAR.pisos-rosubg*((cord(j,lneighbor,2)+cord(j,i,2))/2-PAR.rzbo);
            dlx_l = cord(j,i  ,1)-cord(j,lneighbor,1);
            dly_l = cord(j,i  ,2)-cord(j,lneighbor,2);

            press_norm_r = PAR.pisos-rosubg*((cord(j,rneighbor,2)+cord(j,i,2))/2-PAR.rzbo);
            dlx_r = cord(j,rneighbor,1)-cord(j,i  ,1);
            dly_r = cord(j,rneighbor,2)-cord(j,i  ,2);
            fcx = fcx - 0.5*press_norm_l*dly_l-0.5*press_norm_r*dly_r;
            fcy = fcy + 0.5*press_norm_l*dlx_l+0.5*press_norm_r*dlx_r;

            blx = 1.e+17;
            //write(*,*) i,pisos,force(nz,i,1),force(nz,i,2),press_norm_l,press_norm_r,dlx_l,dlx_r,dly_l,dly_r
        }
    }

    const int ncodx = ncod(j,i,1);
    const int ncody = ncod(j,i,2);

    // BALANCE-OFF
    if( (ncodx & 1) || j<=PAR.n_boff_cutoff )
        blx = 0;
    else
        blx = fabs(fcx) / (blx + 1.e-9);


    if( (ncody & 2) || j<=PAR.n_boff_cutoff )
        bly = 0;
    else
        bly = fabs(fcy) / (bly + 1.e-9);


    // DAMPING
    double vx, vy;
    vx = vel(j,i,1);
    vy = vel(j,i,2);
    if( !(ncodx & 1) && fabs(vx)>1.e-13 ) {
        fcx = fcx - PAR.demf*copysign(fcx, vx);
    }

    if( !(ncody & 2) && fabs(vy)>1.e-13 ) {
        fcy = fcy - PAR.demf*copysign(fcy, vy);
    }

    // VELOCITIES FROM FORCES
    const int iunknown = 0;
    if( ncodx == 1 ) {
        vx = bc(j,i,1) ;
        //            vx = 0.0;

        //        write(*,*) i,j,vx
    }
    else {
        vx = vx + dt*fcx/(amass(j,i)*drat*drat);
    }
    if( ncody == 1 ) {
        vy = bc(j,i,2);
        if(iunknown==1) {
            vy = bc(j,i,2)* sin(time*3.14159/(2*PAR.sec_year));
            //write(*,*) bc(j,i,2), sin(time*3.14159/(2*sec_year));
        }
        //        write(*,*) i,j,vy
    }
    else {
        vy = vy + dt*fcy/(amass(j,i)*drat*drat);
    }

    // Prestress to form the topo when density differences are present WITHOUT PUSHING OR PULLING!
    if (PAR.i_prestress && (time < 200000.*PAR.sec_year)) {
        // node is on bottom/left/right boundary?
        if(j==nz || i==1 || i==nx) {
            vx = 0.0;
            vy = 0.0;
        }
    }

    // Storing the result
    force(j,i,1) = fcx;
    force(j,i,2) = fcy;
    balance(j,i,1) = blx;
    balance(j,i,2) = bly;
    vel(j,i,1) = vx;
    vel(j,i,2) = vy;
    return;
}


extern "C"
void cu_flac(double *force, double *balance, double *vel,
             double *cord, double *stress0, double *temp,
             double *rmass, double *amass,
             const double *bc, const int *ncod,
             const double *time, const double *time_t,
             const double *dtmax_therm, const double *dt,
             const int *nloop, const int *itherm,
             const int *ifreq_rmasses, const int *ifreq_imasses,
             const int *pnx, const int *pnz)
{
    extern void fl_therm_(void);
    extern void fl_srate_(void);
    extern void fl_rheol_(void);
    extern void bc_update_(void);
    extern void fl_move_(void);
    extern void rmasses_(void);
    extern void dt_mass_(void);
    extern void dt_adjust_(void);

    static double *cord_d, *stress0_d, *temp_d, *rmass_d, *force_d, *balance_d, *amass_d, *bc_d, *vel_d;
    static int *ncod_d;
    static int first = 1;

    static cudaStream_t stream1;
    const int nx = *pnx;
    const int nz = *pnz;


    if(first) {
        first = 0;
        //fprintf(stderr, "addr: %d %d %d %d\n", cord, temp, vel, stress0);

        cudaMalloc((void **) &cord_d, nx*nz*2*sizeof(double));
        cudaMalloc((void **) &stress0_d, nx*nz*ntriag*nstr*sizeof(double));
        cudaMalloc((void **) &temp_d, nx*nz*sizeof(double));
        cudaMalloc((void **) &rmass_d, nx*nz*sizeof(double));
        cudaMalloc((void **) &amass_d, nx*nz*sizeof(double));
        cudaMalloc((void **) &bc_d, nx*nz*2*sizeof(double));
        cudaMalloc((void **) &force_d, nx*nz*2*sizeof(double));
        cudaMalloc((void **) &balance_d, nx*nz*2*sizeof(double));
        cudaMalloc((void **) &vel_d, nx*nz*2*sizeof(double));
        cudaMalloc((void **) &ncod_d, nx*nz*2*sizeof(int));

        cudaMemcpyAsync(bc_d, bc, nx*nz*2*sizeof(double),
                        cudaMemcpyHostToDevice, stream1);
        cudaMemcpyAsync(ncod_d, ncod, nx*nz*2*sizeof(int),
                        cudaMemcpyHostToDevice, stream1);


        cudaMemcpyAsync(cord_d, cord, nx*nz*2*sizeof(double),
                        cudaMemcpyHostToDevice, stream1);
        cudaMemcpyAsync(vel_d, vel, nx*nz*2*sizeof(double),
                        cudaMemcpyHostToDevice, stream1);
        cudaMemcpyAsync(rmass_d, rmass, nx*nz*sizeof(double),
                        cudaMemcpyHostToDevice, stream1);
        cudaMemcpyAsync(amass_d, amass, nx*nz*sizeof(double),
                        cudaMemcpyHostToDevice, stream1);

        checkCUDAError("cu_flac: Allocating memory");
    }


    if(*time - *time_t > *dtmax_therm/10) fl_therm_();
    if(*itherm == 2) return;

    cudaMemcpyAsync(temp_d, temp, nx*nz*sizeof(double),
                    cudaMemcpyHostToDevice, stream1);

    fl_srate_();

    fl_rheol_();
    cudaMemcpyAsync(stress0_d, stress0, nx*nz*ntriag*nstr*sizeof(double),
                    cudaMemcpyHostToDevice, stream1);

    if(param.ynstressbc == 1) {
        bc_update_();
        cudaMemcpyAsync(force_d, force, nx*nz*2*sizeof(double),
                        cudaMemcpyHostToDevice, stream1);
        cudaMemcpyAsync(balance_d, balance, nx*nz*2*sizeof(double),
                        cudaMemcpyHostToDevice, stream1);
    }

    dim3 dimBlock(nthx,nthz);
    dim3 dimGrid(nx/nthx+1, nz/nthz+1);
    const double drat = 1.0;//(*dt < *dt_elastic) ? 1.0 : (*dt / *dt_elastic);

    cu_fl_node<<<dimGrid, dimBlock, stream1>>>(force_d, balance_d, vel_d,
                                               cord_d, stress0_d, temp_d,
                                               rmass_d, amass_d,
                                               bc_d, ncod_d,
                                               *dt, drat, *time,
                                               nx, nz);
    cudaStreamSynchronize(stream1);
    cudaMemcpy(vel, vel_d, nx*nz*2*sizeof(double),
               cudaMemcpyDeviceToHost);

    cudaMemcpyAsync(force, force_d, nx*nz*2*sizeof(double),
                    cudaMemcpyDeviceToHost, stream1);
    cudaMemcpyAsync(balance, balance_d, nx*nz*2*sizeof(double),
                    cudaMemcpyDeviceToHost, stream1);

    if(1) {
        fl_move_();
        cudaMemcpyAsync(cord_d, cord, nx*nz*2*sizeof(double),
                        cudaMemcpyHostToDevice, stream1);
    }

    if( (*nloop % *ifreq_rmasses) == 0 ) {
        rmasses_();
        cudaMemcpyAsync(rmass_d, rmass, nx*nz*sizeof(double),
                        cudaMemcpyHostToDevice, stream1);
    }

    if( (*nloop % *ifreq_imasses) == 0 ) {
        dt_mass_();
        cudaMemcpyAsync(amass_d, amass, nx*nz*sizeof(double),
                        cudaMemcpyHostToDevice, stream1);
    }

    dt_adjust_();

    if( (*nloop % *ifreq_rmasses) == 0 )
        cudaMemcpy(rmass, rmass_d, nx*nz*sizeof(double),
                   cudaMemcpyDeviceToHost);

    checkCUDAError("cu_flac: end");

    return;
}


/*
 * Copying Fortran variables to C and CUDA
 */
extern "C" __host__
void cu_copy_param_(int *irheol, double *visc,
                    double *den, double *alfa, double *beta,
                    double *pln, double *acoef, double *eactiv,
                    double *rl, double *rm, double *coha, double *cohdisp,
                    double *phimean, double *phidisp, double *psia,
                    double *conduct, double *cp,
                    double *ts, double *tl, double *tk, double *fk,
                    double *g, double *pisos, double *drosub,
                    double *rzbo, double *demf,
                    double *sec_year, double *ynstressbc,
                    double *dt_scale, double *frac, double *fracm,
                    double *strain_inert, double *vbc,
                    int *lphase,
                    int *nyhydro, int *iphsub,
                    int *n_boff_cutoff, int *i_prestress,
                    int *iint_marker, int *nphasl, int *idt_scale)
{

    memcpy(param.irheol, irheol, 20*sizeof(int));
    memcpy(param.visc, visc, 20*sizeof(double));
    memcpy(param.den, den, 20*sizeof(double));
    memcpy(param.alfa, alfa, 20*sizeof(double));
    memcpy(param.beta, beta, 20*sizeof(double));
    memcpy(param.pln, pln, 20*sizeof(double));
    memcpy(param.acoef, acoef, 20*sizeof(double));
    memcpy(param.eactiv, eactiv, 20*sizeof(double));
    memcpy(param.rl, rl, 20*sizeof(double));
    memcpy(param.rm, rm, 20*sizeof(double));
    memcpy(param.coha, coha, 20*sizeof(double));
    memcpy(param.cohdisp, cohdisp, 20*sizeof(double));
    memcpy(param.phimean, phimean, 20*sizeof(double));
    memcpy(param.phidisp, phidisp, 20*sizeof(double));
    memcpy(param.psia, psia, 20*sizeof(double));
    memcpy(param.conduct, conduct, 20*sizeof(double));
    memcpy(param.cp, cp, 20*sizeof(double));
    memcpy(param.ts, ts, 20*sizeof(double));
    memcpy(param.tl, tl, 20*sizeof(double));
    memcpy(param.tk, tk, 20*sizeof(double));
    memcpy(param.fk, fk, 20*sizeof(double));

    memcpy(param.lphase, lphase, 20*sizeof(int));

    param.g = *g;
    param.pisos = *pisos;
    param.drosub = *drosub;
    param.rzbo = *rzbo;
    param.demf = *demf;
    param.sec_year = *sec_year;
    param.ynstressbc = *ynstressbc;
    param.dt_scale = *dt_scale;
    param.frac = *frac;
    param.fracm = *fracm;
    param.strain_inert = *strain_inert;
    param.vbc = *vbc;


    param.nyhydro = *nyhydro;
    param.iphsub = *iphsub;
    param.n_boff_cutoff = *n_boff_cutoff;
    param.i_prestress = *i_prestress;
    param.iint_marker = *iint_marker;
    param.nphasl = *nphasl;
    param.idt_scale = *idt_scale;

    //fprintf(stderr, "1: %e %e %e %e\n", *demf, *sec_year, *g, *ynstressbc);

    // copy to CUDA constant memory
    cudaMemcpyToSymbol("PAR", &param, sizeof(flac_param_t), 0, cudaMemcpyHostToDevice);
    checkCUDAError("cu_copy_param");
    return;
}
