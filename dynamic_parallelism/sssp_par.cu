#include "sssp_par.cuh"

int main() {
  unsigned int N;
  unsigned int degree;
  unsigned int M;
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
  cudaEvent_t start, stop;
  float elapsedTime;
  cudaEventCreate(&start);
  cudaEventRecord(start,0);
  // execute dijkstra compound frontiers
  DA2CF<<<1,1>>>(c_dev, u_dev, f_dev, e_dev, w_dev, v_dev, N, P);
  cudaEventCreate(&stop);
  cudaEventRecord(stop,0);
  cudaEventSynchronize(stop);

  cudaEventElapsedTime(&elapsedTime, start,stop);
  cout << elapsedTime/1000.0f << " " << N << endl;

  // free allocated memory on the GPU
  cudaFree(c_dev);
  cudaFree(f_dev);
  cudaFree(u_dev);
  cudaFree(v_dev);
  cudaFree(w_dev);
  cudaFree(e_dev);
  return 0;
}
