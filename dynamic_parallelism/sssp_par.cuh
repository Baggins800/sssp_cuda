#include <cuda.h>
#include <iostream>
#include <vector>
#include <ctime>
#include <stdio.h>
#define TPB 512
#define INF 999999999

using namespace std;
struct Node {
  unsigned int start;
  unsigned int adj;
};

__global__ void intialize(unsigned int *c_dev,
                          bool *u_dev,
                          bool *f_dev,
                          unsigned int N, unsigned int start, unsigned int idx) {

  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < N) {
    c_dev[tid + idx] = INF;
    f_dev[tid + idx] = false;
    u_dev[tid + idx] = true;
  }
  __syncthreads();
  if (tid == start) {
    c_dev[tid + idx] = 0;
    f_dev[tid + idx] = true;
    u_dev[tid + idx] = false;
  }
}

__global__ void relax_f(unsigned int *c_dev,
                        bool *u_dev,
                        bool *f_dev,
                        unsigned int *e_dev,
                        unsigned int *w_dev,
                        Node *v_dev,
                        unsigned int N, unsigned int idx) {

  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < N) {
    if (f_dev[tid + idx]) {
      for (int i = v_dev[tid].start; 
           i < v_dev[tid].start + v_dev[tid].adj; 
           i++) {
        unsigned int succ = e_dev[i];
        if (u_dev[succ + idx]) {
          atomicMin(&c_dev[succ + idx], c_dev[tid + idx] + w_dev[i]);
        }
      }
    }
  }
}

__global__ void update(unsigned int * c_dev,
                       bool *f_dev, bool *u_dev,
                       unsigned int mssp, unsigned int N, unsigned int idx) {
  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < N) {
    f_dev[tid + idx] = false;
    if (c_dev[tid + idx] == mssp) {
      u_dev[tid + idx] = false;
      f_dev[tid + idx] = true;
    }
  }
}

__global__ void minimum(unsigned int *c_dev,
                        bool *u_dev, unsigned int &mssp, unsigned int N, unsigned int idx) {
  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < N) {
    if (u_dev[tid + idx] && (c_dev[tid + idx] < mssp)) {
      atomicMin(&mssp, c_dev[tid + idx]);
    }
  }
}

__global__ void DA2CF(unsigned int *c_dev, 
           bool *u_dev, bool *f_dev, 
           unsigned int *e_dev, 
           unsigned int *w_dev,
           Node *v_dev,
           unsigned int N, vector<unsigned int> &P, unsigned int * mssp) {
  unsigned int idx = threadIdx.x + blockDim.x * blockIdx.x;
  if (idx < N) {
    unsigned int extrablock = N % TPB > 0 ? 1 : 0;
    intialize<<<N / TPB + extrablock, TPB>>>(c_dev, u_dev, f_dev, N, idx, idx * N);
    cudaDeviceSynchronize();
    mssp[idx] = 0;
    while (mssp[idx] != INF) {
      mssp[idx] = INF;
      relax_f<<<N / TPB + extrablock, TPB>>>(c_dev, u_dev, f_dev, e_dev, w_dev, v_dev, N, idx * N);
      cudaDeviceSynchronize();
      minimum<<< N / TPB + extrablock, TPB >>>(c_dev, u_dev, mssp[idx], N, idx * N);   
      cudaDeviceSynchronize();
      update<<< N / TPB + extrablock, TPB >>>(c_dev, f_dev, u_dev, mssp[idx], N, idx * N);
      cudaDeviceSynchronize();
    }
  }
}
