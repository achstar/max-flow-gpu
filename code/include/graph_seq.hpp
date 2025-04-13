#include <vector>
#include <string>
using namespace std;

struct Edge
{
    int src;
    int dest;
    int capacity;   // initialize to 0 for reverse edges, change as we go
    int flow = 0; // don't need value for residual flow, it's the same as outgoing flow
};

struct Vertex
{   
    int id;
    int label = 0; 
    int excess = 0; 
    vector<Edge> outgoing_edges;
};

class Graph
{
    vector<Vertex> vertices;
    int source_node;
    int sink_node;
    
    void init_preflow(int s);

    bool push(int u);

    void relabel(int u);

public:

    void addEdge(int src, int dest, int capacity, int flow); // function to add an edge

    int maxFlow(int s, int t); // function that returns maximum flow from source s to sink t
    Graph(int m, int n, int source, int sink);
};
