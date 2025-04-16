#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "graph_seq.hpp"

#define INF 1e9

__global__ void
push_relabel_kernel(int num_nodes, int source, int sink, int* excess, int* labels, int* cf_adj)
{
    unsigned int u = (blockIdx.x * blockDim.x) + threadIdx.x;
    if ((u < num_nodes) && (u != sink))
    {
        int cycle = num_nodes; // replace with KERNEL_CYCLES?
        while (cycle > 0)
        {
            if ((excess[u] > 0) && (labels[u] < num_nodes))
            {
                int e_prime = excess[u];
                int min_label = INF;
                int min_v = 0;
                for (int v = 0; v < num_nodes; v++)
                {
                    if (cf_adj[u*num_nodes + v] > 0)
                    {
                        // printf("Residual at edge (%d, %d) is %d\n", u, v, cf_adj[u*num_nodes + v]);
                        int curr_label = labels[v];
                        if (curr_label < min_label)
                        {
                            min_v = v;
                            min_label = curr_label;
                        }
                    }
                }
                // printf("lowest dest node is %d, with label %d, u is %d\n", min_v, min_label, u);
                if (labels[u] > min_label)
                {
                    // printf("Excess of %d is %d\n", min_v, excess[min_v]);
                    int d = min(e_prime, cf_adj[u*num_nodes + min_v]);
                    atomicAdd(&cf_adj[min_v*num_nodes + u], d);
                    atomicSub(&cf_adj[u*num_nodes + min_v], d);
                    atomicAdd(&excess[min_v], d);
                    atomicSub(&excess[u], d);
                    // printf("Pushing %d at edge (%d, %d)\n", d, u, min_v);
                    // printf("Excess of %d is now updated to %d\n", min_v, excess[min_v]);
                }
                else
                {
                    labels[u] = min_label + 1;
                    // printf("Relabeling node %d to %d\n", u, labels[u]);
                }
            }
            cycle--;
        }
    }
}

void Graph::globalRelabel(int num_nodes, int source, int sink, std::vector<int>& excess, std::vector<int>& labels, std::vector<int>& cf_adj, std::vector<bool>& marked)
{
    // ensure that source label is no more than 1 greater than dest label
    for (int u = 0; u < num_nodes; u++)
    {
        int source_label = labels[u];
        for (const Edge& e : vertices[u].outgoing_edges) {
            int v = e.dest;
            int dest_label = labels[v];
            if (source_label > dest_label + 1)
            {
                int flow = cf_adj[u*num_nodes + v];
                excess[u] -= flow;
                excess[v] += flow;
                cf_adj[v*num_nodes + u] += flow;
                cf_adj[u*num_nodes + v] = 0;
            }
        }
    }

    // BFS from sink to reset labels
    std::vector<bool> visited(N, false);
    std::queue<int> q;

    labels[sink] = 0;
    visited[sink] = true;
    q.push(sink);

    while (!q.empty()) {
        // pop from queue
        int u = q.front(); 
        q.pop();
    
        for (int v = 0; v < num_nodes; v++) {
            if (cf_adj[v*num_nodes + u] > 0 && !visited[v]) { // reverse residual edge exists
                labels[v] = labels[u] + 1; // 1 more than parent node height
                visited[v] = true;
                q.push(v);
            }
        }
    }

    for (int u = 0; u < num_nodes; u++) {
        if (!visited[u] && !marked[u]) { // if not visited and not relabeled
            excess_total -= excess[u];
            marked[u] = true;
            excess[u] = 0;
        }
    }
}

int Graph::maxFlowParallel(int s, int t)
{
    int N = vertices.size();
    int threadsPerBlock = 256;
    int numberOfBlocks = (N / threadsPerBlock) + 1;
    // Host data structures (excess, labels, capacity/flow)
    std::vector<bool> marked(N, false);
    std::vector<int> h_excess(N, 0);
    std::vector<int> h_labels(N, 0);
    // necessary to have "flattened" ds for GPU locality
    std::vector<int> h_cf(N*N, 0); // will need to adjust this for larger test cases because we will run out of memory
    init_preflow(s);
    // printf("Done initializing preflow with excess total %d\n", excess_total);
    for (int u = 0; u < N; u++) {
        h_labels[u] = vertices[u].label;
        h_excess[u] = vertices[u].excess;
        // printf("Node %d has label %d and excess %d\n", u, h_labels[u], h_excess[u]);
        for (const Edge& e : vertices[u].outgoing_edges) {
            h_cf[e.src*N + e.dest] = e.capacity; // "forward" flow
            // No need to add reverse flow since that should already be accounted for in init_preflow?
            // printf("Residual capacity for (%d, %d) is %d\n", e.src, e.dest, e.capacity);
        }
    }

    // Allocate device memory
    int* d_excess, *d_labels, *d_cf;
    cudaMalloc((void **)&d_excess, N*sizeof(int));
    cudaMalloc((void **)&d_labels, N*sizeof(int));
    cudaMalloc((void **)&d_cf, N*N*sizeof(int));

    cudaMemcpy(d_excess, h_excess.data(), N*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_cf, h_cf.data(), N*N*sizeof(int), cudaMemcpyHostToDevice);

    printf("Starting loop\n");
    while (excess_total != h_excess[s] + h_excess[t])
    {
        cudaMemcpy(d_labels, h_labels.data(), N*sizeof(int), cudaMemcpyHostToDevice);
        // printf("Launching kernel\n");
        push_relabel_kernel<<<numberOfBlocks, threadsPerBlock>>>(N, s, t, d_excess, d_labels, d_cf);
        // printf("Launched kernel\n");
        cudaMemcpy(h_excess.data(), d_excess, N*sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_labels.data(), d_labels, N*sizeof(int), cudaMemcpyDeviceToHost);
        // cudaMemcpy(h_cf.data(), d_cf, N*N*sizeof(int), cudaMemcpyDeviceToHost);

        // globalRelabel(N, s, t, h_excess, h_labels, h_cf, marked);
        for(int l : h_labels){
            if(l >= N){
                // printf("setting excess_total now ....\n");
                excess_total = h_excess[s] + h_excess[t];
                break;
            }
        }
        // printf("Excess total: %d\n", excess_total);
        // printf("Excess target: %d\n", h_excess[t]);
    }
    cudaFree(d_excess);
    cudaFree(d_labels);
    cudaFree(d_cf);
    return h_excess[t];
}