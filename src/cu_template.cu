// -*- C++ -*-
#ifndef __CU_TEMPLATE_CU__
#define __CU_TEMPLATE_CU__


#include "cuda.h"


////////////////////////////////////////

enum reduction_op {
    ADD, MAX, MAXABS, MIN, MULTIPLY
};

/*
  Input -- data: array T[M]
  Output-- result: pointer to a scalar
  Must be launched with block=1, blockDim.x is power of 2
 */
template <class T, int M, enum reduction_op op>
__global__
void reduction_kernel(T *data, T *result)
{
    extern __shared__ T smem[];
    unsigned const int tid = threadIdx.x;

    T s = data[tid];
    if(op == MAXABS) s = abs(s);

    for(unsigned int i = tid + blockDim.x; i<M; i+=blockDim.x) {
        switch (op) {
        case ADD:
            s += data[i]; break;
        case MAX:
            s = max(s, data[i]); break;
        case MAXABS:
            s = max(s, abs(data[i])); break;
        case MIN:
            s = min(s, data[i]); break;
        case MULTIPLY:
            s *= data[i]; break;
        }
    }
    smem[tid] = s;
    __syncthreads();

    for(unsigned int j=blockDim.x/2; j>0; j>>=1) {
        if(tid < j) {
            switch (op) {
            case ADD:
                smem[tid] += smem[tid + j]; break;
            case MAX:
                smem[tid] = max(smem[tid], smem[tid + j]); break;
            case MAXABS:
                smem[tid] = max(smem[tid], abs(smem[tid + j])); break;
            case MIN:
                smem[tid] = min(smem[tid], smem[tid + j]); break;
            case MULTIPLY:
                smem[tid] *= smem[tid + j]; break;
            }
        }
        __syncthreads();
    }
    
    // write result for this block to global mem
    if (tid == 0) *result = smem[0];
    return;
}


/*
  Input -- data: array T[M]
 */
template <class T, int M, enum reduction_op op>
__host__
T reduction(T *data_d)
{
    T *res_d, result;
    const int nb = (M > 512) ? 512 : M;
    cudaMalloc((void **) &res_d, sizeof(T));
 
    reduction_kernel<T,M,op><<<1,nb,nb*sizeof(T)>>>(data_d, res_d);
    cudaMemcpy(&result, res_d, sizeof(T), cudaMemcpyDeviceToHost);

    cudaFree(res_d);
    return result;
}


/**** TEST ****/
#if 0
#include <stdio.h>

__host__
void test_reduction()
{
    const int n = 201*501;
    double a[n], *d;
    fprintf(stderr, "1\n");
    
    for(int i=0; i<n; i++) a[i] = -0.1*i;
    cudaMalloc((void **) &d, n*sizeof(double));
    cudaMemcpy(d, a, n*sizeof(double), cudaMemcpyHostToDevice);
    cudaThreadSynchronize();
    fprintf(stderr, "4\n");

    const enum reduction_op op = ADD;
    double result = reduction<double,n,op>(d);
    
    cudaThreadSynchronize();
    fprintf(stderr, "5\n");
    fprintf(stderr, "result = %f\n", result);
    return;
}


int main(int, char**)
{
    test_reduction();
    return 0;
}
#endif

#endif
