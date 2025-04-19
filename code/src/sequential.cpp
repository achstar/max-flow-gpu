#include "graph_seq.hpp"
#include <queue>

// graph constructor to initialize nodes
Graph::Graph(int n, int m, int source, int sink) {
    source_node = source;
    sink_node = sink;
    N = n;
    M = 0; // increment in addEdge
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
    M++;
}

// initializes preflow for edges coming out of source
void Graph::init_preflow(int s)
{
    vertices[s].label++;
    // go through all outgoing edges from source
    for(int i = 0; i < vertices[s].outgoing_edges.size(); i++){
        Edge& curr_edge = vertices[s].outgoing_edges[i];
        Vertex &dest_vertex = vertices[curr_edge.dest];

        // update flow for this edge --> capacity
        curr_edge.flow = curr_edge.capacity;
        curr_edge.capacity = 0;

        // update excess on node we push to
        dest_vertex.excess += curr_edge.flow;

        // take care of reverse edge from dest node
        for(int j = 0; j < dest_vertex.outgoing_edges.size(); j++){
            if(dest_vertex.outgoing_edges[j].dest == s){
                dest_vertex.outgoing_edges[j].capacity += curr_edge.flow;
                excess_total += curr_edge.flow;
            }
        }        
    }
}

// returns whether or not push was successful
bool Graph::push(Vertex& vertex){
    // printf("Current vertex is %d\n", vertex.id);
    for(int i = 0; i < vertex.outgoing_edges.size(); i++){
        Edge& curr_edge = vertex.outgoing_edges[i];
        Vertex& dest_vertex = vertices[curr_edge.dest];
        // printf("{%d, %d} with capacity %d and flow %d\n", curr_edge.src, curr_edge.dest, curr_edge.capacity, curr_edge.flow);
        if(curr_edge.capacity != 0 && vertex.excess > 0){
            // current label is higher, we can push
            if(vertex.label == dest_vertex.label + 1){
                int flow = min(vertex.excess, curr_edge.capacity);
                // printf("Pushing %d flow... from node %d to node %d\n", flow, vertex.id, dest_vertex.id);
                // add flow 
                curr_edge.flow += flow;
                curr_edge.capacity -= flow;
                // subtract excess
                vertex.excess -= flow;
                dest_vertex.excess += flow;
                // update reverse flow 
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
    // printf("Failed to push from node %d\n", vertex.id);
    return false;
}

void Graph::relabel(Vertex& vertex){
    // printf("relabeling vertex %d from %d to %d\n...", vertex.id, vertex.label, vertex.label+1);
    vertex.label++;

}

int Graph::maxFlowSeq(int s, int t)
{
    init_preflow(s);
    printf("done initializing preflow\n");
    // for(int i = 0; i < vertices.size(); i++){
    //     printf("\nVertex %d - label: %d, excess: %d\n", vertices[i].id, vertices[i].label, vertices[i].excess);
    //     for(int j = 0; j < vertices[i].outgoing_edges.size(); j++){
    //         printf("{%d, %d} has cap: %d and flow: %d\n", vertices[i].outgoing_edges[j].src,vertices[i].outgoing_edges[j].dest, vertices[i].outgoing_edges[j].capacity,vertices[i].outgoing_edges[j].flow);
    //     }
    // }

    while(excess_total != vertices[0].excess + vertices[N-1].excess){
        // printf("excess total: %d e(s): %d, e(t): %d\n", excess_total, vertices[0].excess, vertices[N-1].excess );
        for(Vertex& v : vertices){
            if(v.excess != 0 && v.id != t){
                if(v.label > N + 10){
                    excess_total = vertices[0].excess + vertices[N-1].excess;
                    break;
                }
                bool pushed = push(v);
                if (!pushed) {
                    relabel(v);
                }
            }
        }
    }
    return vertices[t].excess;
}
