#include <vector>
#include <string>
#include <queue>
using namespace std;

struct Edge
{
    int src;
    int dest;
    int capacity;   // initialize to 0 for reverse edges, change as we go
    int flow = 0; // don't need value for residual flow, it's the same as outgoing flow
    bool operator==(const Edge& other) const {
        return src == other.src && dest == other.dest;
    }
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
    int N;
    int M;
    vector<Vertex> vertices;
    int source_node;
    int sink_node;
    int excess_total = 0;

    void init_preflow(int s);

    bool push(Vertex& vertex);

    void relabel(Vertex& vertex);
    void globalRelabel(int num_nodes, int source, int sink, std::vector<int>& excess, std::vector<int>& labels, std::vector<int>& cf_adj, std::vector<bool>& marked);

public:

    void addEdge(int src, int dest, int capacity, int flow); // function to add an edge
    int check_excess();
    int maxFlowSeq(int s, int t); // function that returns maximum flow from source s to sink t
    int maxFlowParallel(int s, int t);
    Graph(int n, int m, int source, int sink);
};
