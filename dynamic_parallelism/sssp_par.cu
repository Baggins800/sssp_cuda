#include "sssp_par.cuh"

int main() {
  unsigned int N;
  unsigned int degree;
  unsigned int M;
  cin >> N;
  cin >> degree;
  cin >> M;
  vector<Node> v_host(N);
  vector<unsigned int> e_host(M);
  vector<unsigned int> w_host(M);
  vector<unsigned int> P;
  unsigned int *c_dev, *w_dev, *e_dev, *c_host, * mssp;
  bool *f_dev, *u_dev;
  Node *v_dev;
  c_host = (unsigned int*)malloc(N * N * sizeof(unsigned int));
  vector<int> freq_table(N);
  for (unsigned int i = 0; i < M; i++) {
    unsigned int ia, ib, w;
    cin >> ia >> ib >> w;
    e_host[i] = ib;
    w_host[i] = w;
    freq_table[ia]++;
  }
  v_host[0].start = 0;
  v_host[0].adj = freq_table[0];

  for (unsigned int i = 1; i < N; i++) {
     v_host[i].start = v_host[i - 1].start + v_host[i - 1].adj;
     v_host[i].adj = freq_table[i];
  }
  // allocate frontiers, unresolved and cost vectors on the GPU

  cudaMalloc( (void**)&c_dev, N * N * sizeof(unsigned int) ); 
  cudaMalloc( (void**)&f_dev, N * N * sizeof(bool) ); 
  cudaMalloc( (void**)&u_dev, N * N * sizeof(bool) );
  cudaMalloc( (void**)&v_dev, N * sizeof(Node) );
  cudaMalloc( (void**)&e_dev, M * sizeof(unsigned int) );
  cudaMalloc( (void**)&w_dev, M * sizeof(unsigned int) );
  cudaMalloc( (void**)&mssp, N * sizeof(unsigned int) );

  // copy data to GPU memory
  cudaMemcpy( v_dev, v_host.data(), N * sizeof(Node), cudaMemcpyHostToDevice);
  cudaMemcpy( e_dev, e_host.data(), M * sizeof(unsigned int), cudaMemcpyHostToDevice);
  cudaMemcpy( w_dev, w_host.data(), M * sizeof(unsigned int), cudaMemcpyHostToDevice);
  cudaEvent_t start, stop;
  float elapsedTime;
  cudaEventCreate(&start);
  cudaEventRecord(start,0);
  // execute dijkstra compound frontiers
  int extra_block = N % TPB > 0 ? 1 : 0;
  DA2CF<<<N / TPB + extra_block, TPB>>>(c_dev, u_dev, f_dev, e_dev, w_dev, v_dev, N, P, mssp);
  cudaEventCreate(&stop);
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start,stop);
  cout << elapsedTime/1000.0f << " " << N << endl;
  cudaMemcpy( c_host, c_dev, N * N * sizeof(unsigned int), cudaMemcpyDeviceToHost);
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      cout << c_host[i + N * j] << "\t";
    }
    cout << endl;
  }
  // free allocated memory on the GPU
  cudaFree(c_dev);
  cudaFree(f_dev);
  cudaFree(u_dev);
  cudaFree(v_dev);
  cudaFree(w_dev);
  cudaFree(e_dev);
  cudaFree(mssp);
  free(c_host);
  return 0;
}
