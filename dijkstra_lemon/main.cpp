#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/dijkstra.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
using namespace std;
using namespace lemon;
void read_graph(std::string filename);
typedef ListGraph Graph;
typedef Graph::Node Node;
typedef Graph::Edge Edge;
typedef Graph::EdgeMap<int> LengthMap;
Graph g;
LengthMap len(g);
vector<Node> nodes;
vector<Edge> edges;
void read_() {
  int num_nodes;
  int degree;
  int num_edges;
  cin >> num_nodes;
  cin >> degree;
  cin >> num_edges;
  for (int z = 0; z < num_nodes; z++) {
    nodes.push_back(g.addNode());
  }
  for (int i = 0; i < num_edges; i++) {
    int indexA, indexB, weight;
    cin >> indexA >> indexB >> weight;
    edges.push_back(g.addEdge(nodes[indexA], nodes[indexB]));
    len.set(edges[i], weight);

  }
}
int main (int, char*[]) {
  int begin = 0;
  int end = 0;
  read_();
  end = nodes.size() - 1;
  
  clock_t begintime = clock();
  Dijkstra<Graph, LengthMap> dijkstra_test(g,len);

  dijkstra_test.run(nodes[begin]);
  clock_t endtime = clock();
  double elapsed_time = (double)(endtime - begintime) / CLOCKS_PER_SEC;
  std::cout << elapsed_time  << " " << end + 1 << endl;//<< " " << dijkstra_test.dist(nodes[end]) << std::endl;

/*for (Node v=nodes[end];v != nodes[begin]; v=dijkstra_test.predNode(v)) {
  std::cout << g.id(v) << "<-";
}*/
  return 0;
}
