#include "graph_seq.hpp"

Graph::Graph(int n, int m, int source, int sink) {
    source_node = source;
    sink_node = sink;
    for (int i = 0; i < n; i++) {
        vertices.push_back({i});
    }
}

void Graph::addEdge(int src, int dest, int capacity, int flow)
{
    Edge new_edge = {src, dest, capacity, flow};
    edges.push_back(new_edge);
}

// initializes preflow for edges coming out of source
void Graph::init_preflow(int s)
{
    for(int i = 0; i < edges.size(); i++){
        // edge comes out of source 
        Edge& curr_edge = edges[i];
        if(curr_edge.src == s){
            curr_edge.flow = curr_edge.capacity;
            vertices[curr_edge.dest].excess += curr_edge.flow;
            // add reverse edge? 
            addEdge(curr_edge.dest, curr_edge.src, -curr_edge.flow);
        }
    }
}
// returns whether or not push was successful
bool Graph::push(Vertex vertex){
    for(int i = 0; i < vertex.outgoing_edges.size(); i++){
        Edge& curr_edge = vertex.outgoing_edges[i];
        if(curr_edge.flow < curr_edge.capacity){
            // current label is higher, we can push
            if(vertex.label > curr_edge.dest.label){
                // min(excess, cf)
                int flow = min(vertex.excess, curr_edge.capacity - curr_edge.flow);
                // add flow 
                curr_edge.flow += flow;
                // subtract excess
                vertex.excess -= flow;
                curr_edge.dest.excess += flow;

                // update reverse flow 
                



            }
        }
    }
    

}

void relabel(int u);

void updateResidual(int i, int flow);


