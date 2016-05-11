#include <iostream>
#include <time.h>
#include <cstdlib>
using namespace std;
int main() {

  unsigned int nodes, degree;
  cin >> nodes;
  cin >> degree;
  unsigned int sumi = 0;
  for (int z = 1; z < degree; z++) {
    sumi += z;
  }

  srand(time(NULL));
  unsigned int num_e = nodes * degree - degree * degree + sumi;
  cout << nodes << endl;
  cout << degree << endl;
  cout << num_e << endl;
  for (int i = 0; i < nodes; i++) {
    for (int a = 1; a < degree + 1; a++) {
      if (a + i < nodes) {
        cout << i << " " << i + a << " " << rand() % 10 + 1<< endl;
      }
    }
  }
  return 0;
}
