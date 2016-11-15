#include <cuda.h>
#include <iostream>
#include <vector>
#include <ctime>
#include <stdio.h>
#define TPB 32
#define INF 99999999

using namespace std;
struct Node {
  unsigned int start;
  unsigned int adj;
};

__global__ void intialize(unsigned int *c_dev,
                          bool *u_dev,
                          bool *f_dev,
                          unsigned int N) {

  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < N && tid > 0) {
    c_dev[tid] = INF;
    f_dev[tid] = false;
    u_dev[tid] = true;
  }
  if (tid == 0) {
    c_dev[tid] = 0;
    f_dev[tid] = true;
    u_dev[tid] = false;
  }
}

__global__ void relax_adj(unsigned int *c_dev, 
                          Node * v_dev, 
                          unsigned int *e_dev,
                          unsigned int *w_dev,
                          bool *u_dev,
                          unsigned int tid, 
                          unsigned int istart, 
                          unsigned int N) {
  unsigned int i = blockIdx.x * blockDim.x + threadIdx.x;
  if ( i < v_dev[tid].adj) {
    unsigned int succ = e_dev[i + istart];
    if (u_dev[succ]) {
      atomicMin(&c_dev[succ], c_dev[tid] + w_dev[i + istart]);
    }
  }
 
}
__global__ void relax_f(unsigned int *c_dev,
                        bool *u_dev,
                        bool *f_dev,
                        unsigned int *e_dev,
                        unsigned int *w_dev,
                        Node *v_dev,
                        unsigned int N) {

  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < N) {
    if (f_dev[tid]) {
      for (int i = v_dev[tid].start; 
           i < v_dev[tid].start + v_dev[tid].adj; 
           i++) {
        unsigned int succ = e_dev[i];
        if (u_dev[succ]) {
          atomicMin(&c_dev[succ], c_dev[tid] + w_dev[i]);
        }
      }
    }
  }
}

__global__ void update(unsigned int * c_dev,
                       bool *f_dev, bool *u_dev,
                       unsigned int mssp, unsigned int N) {
  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < N) {
    f_dev[tid] = false;
    if (c_dev[tid] == mssp) {
      u_dev[tid] = false;
      f_dev[tid] = true;
    }
  }
}

__global__ void minimum(unsigned int *c_dev,
                        bool *u_dev, unsigned int &mssp, unsigned int N) {
  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < N) {
    if (u_dev[tid] && (c_dev[tid] < mssp)) {
      atomicMin(&mssp, c_dev[tid]);
    }
  }
}

__device__ unsigned int mssp = 0;
__global__ void DA2CF(unsigned int *c_dev, 
           bool *u_dev, bool *f_dev, 
           unsigned int *e_dev, 
           unsigned int *w_dev,
           Node *v_dev,
           unsigned int N, vector<unsigned int> &P) {
  unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;
  if (idx < N) {
    unsigned int extrablock = N % TPB > 0 ? 1 : 0;
    intialize<<<N / TPB + extrablock, TPB>>>(c_dev + idx * N, u_dev, f_dev, N);
    cudaDeviceSynchronize();
    while (mssp != INF) {
      mssp = INF;
      relax_f<<<N / TPB + extrablock, TPB>>>(c_dev + idx * N, u_dev, f_dev, e_dev, w_dev, v_dev, N);
      cudaDeviceSynchronize();
      minimum<<< N / TPB + extrablock, TPB >>>(c_dev + idx * N, u_dev, mssp, N);   
      cudaDeviceSynchronize();
      update<<< N / TPB + extrablock, TPB >>>(c_dev + idx * N, f_dev, u_dev, mssp, N);
      cudaDeviceSynchronize();
    }
  }
}
