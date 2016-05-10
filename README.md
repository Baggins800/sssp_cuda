**SSSP implementation with CUDA**
=================================
Hardware Requirements
---------------------
```
Nvidia device with cuda compute capability of atleast 3.5 (Kepler, Maxwell and better).
CPU with x64 architecture.
Atleast 4096MiB host Random access memory.
Atleast 20GiB free harddrive space.
```
Software Requirements
---------------------
Linux is required to use the scripts.
The experiments were done with Arch Linux.
```
sudo pacman -S evince
sudo pacman -S texlive-most
sudo pacman -S nvidia nvidia-utils cuda
```
Also, install lemon by following the instructions on http://lemon.cs.elte.hu/trac/lemon.
Background
----------
In mathematics and computer science, graph theory is the study of graphs, which are mathematical structures used to model pairwise relations between objects. A graph in this context is made up of vertices, nodes, or points which are connected by edges, arcs, or lines. A graph may be undirected, meaning that there is no distinction between the two vertices associated with each edge, or its edges may be directed from one vertex to another; see Graph (discrete mathematics) for more detailed definitions and for other variations in the types of graph that are commonly considered. Graphs are one of the prime objects of study in discrete mathematics.

In graph theory, the shortest path problem is the problem of finding a path between two vertices (or nodes) in a graph such that the sum of the weights of its constituent edges is minimized.

Dijkstra's algorithm is an algorithm for finding the shortest paths between nodes in a graph, which may represent, for example, road networks. It was conceived by computer scientist Edsger W. Dijkstra in 1956 and published three years later.
Input
-----------------
Run the `experiment` script with parameters `MAX`, `DEGREE`, `STEPS`.
The command `sh experiment 10 3 1000` will allow a graph of 1000; 2000; ...; 10000 nodes and a degree of 3 arcs per node to be generated and processed.
* The source node of the graph is taken as the first node which is index `0`.
* The destination node of the graph is the last node of the graph which is index `number of nodes minus one`
* The time the algorithm takes, with the number of nodes in the graph is stored in the `results` directory.

### Output format for the files

```
P $Q_n$
0.1 1000
0.2 2000
...
1.0 10000
```

Results will be generated and displayed with evince.
In the main directory, the graphs are stored in main.pdf
