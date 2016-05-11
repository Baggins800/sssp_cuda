#include <iostream>
#include <lemon/list_graph.h>
#include <lemon/dijkstra.h>
#include <string.h>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include<lemon/time_measure.h>
#include <chrono>

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
  //typedef std::chrono::high_resolution_clock clock;

  std::chrono::time_point<std::chrono::system_clock> start1, end1;
  start1 = std::chrono::high_resolution_clock::now();

  Dijkstra<Graph, LengthMap> dijkstra_test(g,len);
  dijkstra_test.run(nodes[begin]);

  end1 = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed_seconds = end1-start1;
  //std::time_t end_time = std::chrono::system_clock::to_time_t(end1);

  std::cout << elapsed_seconds.count() << " " << end + 1 << endl;//<< " " << dijkstra_test.dist(nodes[end]) << std::endl;
/*for (Node v=nodes[end];v != nodes[begin]; v=dijkstra_test.predNode(v)) {
  std::cout << g.id(v) << "<-";
}*/
  return 0;
}
