#include "graph_seq.hpp"
#include <queue>
Graph::Graph(int n, int source, int sink) {
    source_node = source;
    sink_node = sink;
    for (int i = 0; i < n; i++) {
        Vertex v;
        v.id = i;
        vertices.push_back(v);
    }
}

void Graph::addEdge(int src, int dest, int capacity, int flow)
{
    Edge new_edge = {src, dest, capacity, flow};
    vertices[src].outgoing_edges.push_back(new_edge);
}

// initializes preflow for edges coming out of source
void Graph::init_preflow(int s)
{
    // go through all outgoing edges from source

    for(int i = 0; i < vertices[s].outgoing_edges.size(); i++){
        Edge& curr_edge = vertices[s].outgoing_edges[i];
        // update flow for this edge
        curr_edge.flow = curr_edge.capacity;

        Vertex &dest_vertex = vertices[curr_edge.dest];
        // update excess on node we push to
        dest_vertex.excess += curr_edge.flow;

        // take care of reverse edge from dest node
        for(int j = 0; j < dest_vertex.outgoing_edges.size(); j++){
            if(dest_vertex.outgoing_edges[j].dest == s){
                dest_vertex.outgoing_edges[j].capacity += curr_edge.flow;
            }
        }        
    }
}

// returns whether or not push was successful
bool Graph::push(Vertex& vertex){
    for(int i = 0; i < vertex.outgoing_edges.size(); i++){
        Edge& curr_edge = vertex.outgoing_edges[i];
        Vertex& dest_vertex = vertices[curr_edge.dest];
        if(curr_edge.flow < curr_edge.capacity){
            // current label is higher, we can push
            if(vertex.label > dest_vertex.label){
                // min(excess, cf)
                int flow = min(vertex.excess, curr_edge.capacity - curr_edge.flow);
                // add flow 
                curr_edge.flow += flow;
                curr_edge.capacity -= flow;
                // subtract excess
                vertex.excess -= flow;
                dest_vertex.excess += flow;

                // update reverse flow 
                // get corresponding edge
                for(int j = 0; j < dest_vertex.outgoing_edges.size(); j++){
                    if(dest_vertex.outgoing_edges[j].dest == vertex){
                        // add to capacity of residual graph edge
                        dest_vertex.outgoing_edges[j].capacity += flow;
                    }
                }
            }
        }
    } 
}

void Graph::relabel(Vertex& vertex){
    // find adjacent node with lowest height + we can push to
    int min_height = N + 5;
    for(int i = 0; i < vertex.outgoing_edges.size(); i++)
    {
        if(vertex.outgoing_edges[i].flow < vertex.outgoing_edges[i].capacity){
            Vertex& dest_vertex = vertex.outgoing_edges[i].dest;
            if(dest_vertex.height < min_height)
            {
                min_height = dest_vertex.height;
            }
        }
    }
    vertex.height = min_height + 1; 

}

int Graph::maxFlow(int s, int t)
{
    bool excess = true;
    init_preflow(s);
    queue<int> active;

    for(Vertex& v : vertices) {
        if(v.id != s && v.id != t) {
            active.push(v.id);
        }
    }

    while(!active.empty()) {
        Vertex& v  = vertices[active.front()];
        bool pushed = push(v);

        if (!pushed) {
            relabel(v);
        }

        if (v.excess == 0) {
            active.pop();
        } else {
            // Move to back for fairness
            active.pop();
            active.push(v);
        }
    }

    return vertices[t].excess;
}
