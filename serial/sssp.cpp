#include <iostream>
#include <stdio.h>
#include <vector>
#include <ctime>
#define inf 999999
#define THREADS_PER_BLOCK 1024
using namespace std;
struct Node {
  unsigned int start;
  unsigned int num;
  int pred;
};
void initialise_serial(vector<unsigned int> &c, vector<bool>& u) {
  for (int z = 0; z < c.size(); z++) {
    c[z] = inf;
    u[z] = true;
  }
  c[0] = 0;
  u[0] = false;
}

void relax_serial(vector<unsigned int>& c, 
                  vector<bool> u, 
                  vector<unsigned int> w,
                  vector<unsigned int> e,
                  unsigned int f, vector<Node>& n) {
  for (int i = 0; i < u.size(); i++) {
    if (u[i]) {
      unsigned int w0 = inf;
      for (int z = n[f].start; z < n[f].start + n[f].num; z++) {
        if (i == e[z]) {
          w0 = w[z];
        }
      }
      if (c[i] > c[f] + w0)
        n[i].pred = f;
      c[i] = min(c[i], c[f] + w0); 
    }
  }
}
void update_serial(vector<unsigned int> c, vector<bool> &u, unsigned int &f,
    unsigned int &mssp) {
  for (int i = 0; i < u.size(); i++) {
    if (u[i]) {
      if (c[i] < mssp) {
        mssp = c[i];
        f = i;
      }
    }
  }
  u[f] = false;
}
int main() {
  unsigned int a, b, degree;
  cin >> a;
  cin >> degree;
  cin >> b;
  vector<unsigned int> e_host(b);
  vector<unsigned int> w_host(b);
  vector<Node> v_host(a);
  vector<unsigned int> c_host(a);
  vector<bool> u_host(a);
  for (int o = 0; o < a; o++) {
    v_host[o].num = 0;
    v_host[o].start = o * degree;
    v_host[o].pred = -1;
  }
  
  for (int i = 0; i < b; i++) {
    int s, d, w;
    cin >> s >> d >> w;
    w_host[i] = w;
    e_host[i] = d;
    v_host[s].num++;
  }

  clock_t begin = clock();
  unsigned int f = 0;
  unsigned int mssp = 0;
  initialise_serial(c_host, u_host);
  while (mssp != inf) {
    relax_serial(c_host, u_host, w_host, e_host, f, v_host);
    mssp = inf;
    update_serial(c_host, u_host, f, mssp);
  }
  clock_t end = clock();
  double time = (double)(end - begin) / CLOCKS_PER_SEC;
  /*int notstart = a - 1;
  int val = 0;
  cout << a - 1 << "<-";
  while ( notstart != 0) {
    cout << v_host[notstart].pred << "<-" << endl;

    for (int i = v_host[v_host[notstart].pred ].start; i < v_host[v_host[notstart].pred ].start +
        v_host[v_host[notstart].pred ].num; i++) {
      if (e_host[i] == notstart)
        val += w_host[i];

          cout << e_host[i] << " "  << w_host[i] << endl;
    }
    notstart = v_host[notstart].pred;
  }*/
 // cout << endl;
  printf("%f %d\n", time, c_host.size());
  return 0;
}
