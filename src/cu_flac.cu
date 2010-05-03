// -*- C++ -*-

#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>

#include "cu_flac.h"

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
