#include <cuda.h>
#include <iostream>
#include <vector>
#include <ctime>
#include <stdio.h>
#include <thrust/host_vector.h>
#include <thrust/extrema.h>
#include <thrust/device_ptr.h>
#include <thrust/device_malloc.h>
#include <thrust/device_free.h>
#include <thrust/device_vector.h>
#define TPB 1024
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
                       unsigned int *mssp, unsigned int N) {
  unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < N) {
    f_dev[tid] = false;
    if (c_dev[tid] == mssp[0]) {
      u_dev[tid] = false;
      f_dev[tid] = true;
    }
  }
}

__global__ void minimum(unsigned int *c_dev,
                        bool *u_dev, 
                        unsigned int *mssp, 
                        unsigned int N) {
  __shared__ unsigned int sdata[TPB];
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x * 2 * blockDim.x + threadIdx.x;
  unsigned int j = i + blockDim.x;
  sdata[tid] = INF;
  if ((j < N) && (i < N)) {
    sdata[tid] = INF;
    unsigned int A = u_dev[i] ? c_dev[i] : INF;
    unsigned int B = u_dev[j] ? c_dev[j] : INF;
    sdata[tid] = min(A, B);
    __syncthreads();
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
      if (tid < s) {
        sdata[tid] = min(sdata[tid], sdata[tid + s]);
      }
      __syncthreads();
    }
    if (tid == 0) {
      mssp[blockIdx.x] = sdata[0];
    }
  }
}

void DA2CF(unsigned int *c_dev, 
           bool *u_dev, bool *f_dev, 
           unsigned int *e_dev, 
           unsigned int *w_dev,
           Node *v_dev,
           unsigned int N, vector<unsigned int> &P) {
  unsigned int extrablock = N % TPB > 0 ? 1 : 0;
  unsigned int blocks = (extrablock + N / TPB) / 2;
  unsigned int *mssp_dev;
  unsigned int *dev_min_list;
  cudaMalloc( (void**)&mssp_dev, sizeof(unsigned int) );
  cudaMalloc( (void**)&dev_min_list, blocks * sizeof(unsigned int) );
  intialize<<<N / TPB + extrablock, TPB>>>(c_dev, u_dev, f_dev, N);
  unsigned int mssp = 0;
  while (mssp != INF) {
    mssp = INF;
    relax_f<<< N / TPB + extrablock, TPB >>>(c_dev, u_dev, f_dev, e_dev, w_dev, v_dev, N);
    minimum<<< blocks, TPB >>>(c_dev, u_dev, dev_min_list, N);
    thrust::device_ptr<unsigned int> dev_ptr = thrust::device_pointer_cast(dev_min_list);
    thrust::device_ptr<unsigned int> min_ptr = 
    thrust::min_element(dev_ptr, dev_ptr + blocks);
    unsigned int min_value[1];
    min_value[0] = min_ptr[0];
    mssp = min_value[0];
    cudaMemcpy( mssp_dev, &min_value, sizeof(unsigned int), cudaMemcpyHostToDevice);
    update<<< N / TPB + extrablock, TPB >>>(c_dev, f_dev, u_dev, mssp_dev, N);
  }
  cudaFree(&dev_min_list);
  cudaFree(&mssp_dev);
}

int main() {
  unsigned int N, degree, M;
  cin >> N;
  cin >> degree;
  cin >> M;
  vector<unsigned int> c_host(N);
  vector<Node> v_host(N);
  vector<unsigned int> e_host(M);
  vector<unsigned int> w_host(M);
  vector<unsigned int> P;
  unsigned int *c_dev, *w_dev, *e_dev;
  bool *f_dev, *u_dev;
  Node *v_dev;
  for (unsigned int i = 0; i < M; i++) {
    unsigned int ia, ib, w;
    cin >> ia >> ib >> w;
    e_host[i] = ib;
    w_host[i] = w;
  }
  for (unsigned int i = 0; i < N - degree; i++) {
     v_host[i].start = i * degree;
     v_host[i].adj = degree;
  }
  unsigned int dec = degree - 1;
  for (unsigned int z = N - degree; z < N; z++) {
    v_host[z].start = v_host[z - 1].start + dec + 1;
    v_host[z].adj = dec;
    dec--;
  }

  // allocate frontiers, unresolved and cost vectors on the GPU
  cudaMalloc( (void**)&c_dev, N * sizeof(unsigned int) ); 
  cudaMalloc( (void**)&f_dev, N * sizeof(bool) ); 
  cudaMalloc( (void**)&u_dev, N * sizeof(bool) );
  cudaMalloc( (void**)&v_dev, N * sizeof(Node) );
  cudaMalloc( (void**)&e_dev, M * sizeof(unsigned int) );
  cudaMalloc( (void**)&w_dev, M * sizeof(unsigned int) );

  // copy data to GPU memory
  cudaMemcpy( v_dev, v_host.data(), N * sizeof(Node), cudaMemcpyHostToDevice);
  cudaMemcpy( e_dev, e_host.data(), M * sizeof(unsigned int), cudaMemcpyHostToDevice);
  cudaMemcpy( w_dev, w_host.data(), M * sizeof(unsigned int), cudaMemcpyHostToDevice);

  // execute dijkstra compound frontiers
  cudaEvent_t start, stop;
  float elapsedTime;
  cudaEventCreate(&start);
  cudaEventRecord(start, 0);
  DA2CF(c_dev, u_dev, f_dev, e_dev, w_dev, v_dev, N, P);
  cudaEventCreate(&stop);
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  cout << elapsedTime / 1000.0f << " " << N << endl;

  // free allocated memory on the GPU
  cudaFree(c_dev);
  cudaFree(f_dev);
  cudaFree(u_dev);
  cudaFree(v_dev);
  cudaFree(w_dev);
  cudaFree(e_dev);
  return 0;
}
