#include <iostream>
#include <time.h>
#include <cstdlib>
#include <vector>
#include <tuple>
#include <algorithm>
using namespace std;
bool sort_first(tuple<int, int, int> a, tuple<int, int, int> b) {
  if (get<0>(a) < get<0>(b)) return true;
  return false; 
}
int main() {
  unsigned int nodes, degree, undirected;
  cin >> nodes;
  cin >> degree;
  cin >> undirected;
  unsigned int sumi = 0;
  for (int z = 1; z < degree; z++) {
    sumi += z;
  }
  srand(time(NULL));
  unsigned int num_e = 2 * (nodes * degree - degree * degree + sumi);
  if (undirected == 0) num_e /= 2;
  cout << nodes << endl;
  cout << degree << endl;
  cout << num_e << endl;
  vector<tuple<int, int, int>> val; 
  for (int i = 0; i < nodes; i++) {
    for (int a = 1; a < degree + 1; a++) {
      if (a + i < nodes) {
        int w = rand() % 10 + 1;
        val.emplace_back(i, i + a, w);
        if (undirected == 1) 
          val.emplace_back(i + a, i, w);
      }
    }
  }
  sort(val.begin(), val.end(), sort_first); 
  for (auto t : val) {
    cout << get<0>(t) << " " << get<1>(t) << " " <<  get<2>(t)<< endl;
  }
  return 0;
}
