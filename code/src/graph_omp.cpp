#include <omp.h>
#include <graph_seq.hpp>
int Graph::maxFlowOmp(int s, int t)
{
    init_preflow(s);
    printf("done initializing preflow\n");
    while (excess_total != vertices[s].excess + vertices[t].excess) {
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < vertices.size(); i++) {
            Vertex& v = vertices[i];
            if (v.excess != 0 && v.id != t) {
                bool pushed = false;
                for (Edge& curr_edge : v.outgoing_edges) {
                    Vertex& dest_vertex = vertices[curr_edge.dest];
                    if (curr_edge.capacity != 0 && v.excess > 0) {
                        if(v.label == dest_vertex.label + 1){
                            int flow = std::min(v.excess, curr_edge.capacity);

                            #pragma omp critical
                            {   
                                curr_edge.flow += flow;
                                curr_edge.capacity -= flow;
                                v.excess -= flow;
                                dest_vertex.excess += flow;

                            // Update reverse edge
                                for (Edge& rev_edge : dest_vertex.outgoing_edges) {
                                    if (rev_edge.dest == v.id) {
                                        rev_edge.capacity += flow;
                                        break;
                                    }
                                }   
                            }                            
                            pushed = true;
                        }
                    }
                }

                if (!pushed) {
                    if(v.id == s){
                        excess_total = vertices[s].excess + vertices[t].excess;
                    }
                    else{
                        #pragma omp atomic
                        v.label++;
                    }
                }
            }
        }
    }

    return vertices[t].excess;
}