#include "graph_seq.hpp"
#include <queue>
Graph::Graph(int n, int m, int source, int sink) {
    source_node = source;
    sink_node = sink;
    for (int i = 1; i < n + 1; i++) {
        Vertex v;
        v.id = i;
        vertices.push_back(v);
    }
}

void Graph::addEdge(int src, int dest, int capacity, int flow)
{
    Edge new_edge = {src, dest, capacity, flow};
    vertices[src-1].outgoing_edges.push_back(new_edge);
}

// initializes preflow for edges coming out of source
void Graph::init_preflow(int s)
{
    // go through all outgoing edges from source
    int idx = s - 1;
    for(int i = 0; i < vertices[idx].outgoing_edges.size(); i++){
        Edge& curr_edge = vertices[idx].outgoing_edges[i];
        // update flow for this edge
        curr_edge.flow = curr_edge.capacity;

        Vertex &dest_vertex = vertices[curr_edge.dest-1];
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
    printf("Current vertex is %d\n", vertex.id);
    for(int i = 0; i < vertex.outgoing_edges.size(); i++){
        Edge& curr_edge = vertex.outgoing_edges[i];
        Vertex& dest_vertex = vertices[curr_edge.dest-1];
        printf("{%d, %d} with capacity %d and flow %d\n", curr_edge.src, curr_edge.dest, curr_edge.capacity, curr_edge.flow);
        if(curr_edge.flow < curr_edge.capacity){
            // current label is higher, we can push
            if(vertex.label > dest_vertex.label){
                // min(excess, cf)
                int flow = min(vertex.excess, curr_edge.capacity - curr_edge.flow);
                printf("Pushing %d flow...\n", flow);
                // add flow 
                curr_edge.flow += flow;
                curr_edge.capacity -= flow;
                // subtract excess
                vertex.excess -= flow;
                dest_vertex.excess += flow;
                printf("dest: %d", dest_vertex.id);
                // update reverse flow 
                // get corresponding edge
                for(int j = 0; j < dest_vertex.outgoing_edges.size(); j++){
                    if(dest_vertex.outgoing_edges[j].dest == vertex.id){
                        // add to capacity of residual graph edge
                        dest_vertex.outgoing_edges[j].capacity += flow;
                    }
                }
                return true;
            }
        }
    } 
    return false;
}

void Graph::relabel(Vertex& vertex){
    // find adjacent node with lowest height + we can push to
    printf("Relabeling vertex %d...\n", vertex.id);
    int min_height = N + 5;
    for(int i = 0; i < vertex.outgoing_edges.size(); i++)
    {
        if(vertex.outgoing_edges[i].flow < vertex.outgoing_edges[i].capacity){
            Vertex& dest_vertex = vertices[vertex.outgoing_edges[i].dest-1];
            if(dest_vertex.label < min_height)
            {
                min_height = dest_vertex.label;
            }
        }
    }
    if(min_height != N + 5){
        vertex.label++; 
    }
    printf("Relabeled vertex %d to %d \n", vertex.id, vertex.label);


}

int Graph::maxFlow(int s, int t)
{
    // for(int i = 0; i < vertices.size(); i++){
    //     printf("\nVertex %d - label: %d, excess: %d\n", vertices[i].id, vertices[i].label, vertices[i].excess);
    //     for(int j = 0; j < vertices[i].outgoing_edges.size(); j++){
    //         printf("{%d, %d} has cap: %d and flow: %d\n", vertices[i].outgoing_edges[j].src,vertices[i].outgoing_edges[j].dest, vertices[i].outgoing_edges[j].capacity,vertices[i].outgoing_edges[j].flow);
    //     }
    // }
    init_preflow(s);
    printf("done initializing preflow\n");

    queue<int> active;

    for(Vertex& v : vertices) {
        if(v.id != s && v.id != t) {
            active.push(v.id);
        }
    }
    printf("done initializing excess queue\n");

    while(!active.empty()) {
        Vertex& v  = vertices[active.front()-1];
        printf("pushing from vertex %d\n", v.id);
        bool pushed = push(v);
        printf("done push\n");
        for(int i = 0; i < vertices.size(); i++){
            printf("\nVertex %d - label: %d, excess: %d\n", vertices[i].id, vertices[i].label, vertices[i].excess);
            for(int j = 0; j < vertices[i].outgoing_edges.size(); j++){
                printf("{%d, %d} has cap: %d and flow: %d\n", vertices[i].outgoing_edges[j].src,vertices[i].outgoing_edges[j].dest, vertices[i].outgoing_edges[j].capacity,vertices[i].outgoing_edges[j].flow);
            }
        }


        if (!pushed) {
            relabel(v);
            for(int i = 0; i < vertices.size(); i++){
                printf("\nVertex %d - label: %d, excess: %d\n", vertices[i].id, vertices[i].label, vertices[i].excess);
                for(int j = 0; j < vertices[i].outgoing_edges.size(); j++){
                    printf("{%d, %d} has cap: %d and flow: %d\n", vertices[i].outgoing_edges[j].src,vertices[i].outgoing_edges[j].dest, vertices[i].outgoing_edges[j].capacity,vertices[i].outgoing_edges[j].flow);
                }
            }
        }

        if (v.excess == 0) {
            active.pop();
        } else {
            // Move to back for fairness
            active.pop();
            active.push(v.id);
        }
    }

    return vertices[t-1].excess;
}
