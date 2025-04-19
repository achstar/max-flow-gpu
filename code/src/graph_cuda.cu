#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "graph_seq.hpp"

#define INF 1e9

__global__ void
push_relabel_kernel(int num_nodes, int source, int sink, int* excess, int* labels, int* cf,  int* edge_starts, int* edge_dests, int* reverse_edge_index)
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
                // for (int v = 0; v < num_nodes; v++)
                // {
                //     if (cf_adj[u*num_nodes + v] > 0)
                //     {
                //         // printf("Residual at edge (%d, %d) is %d\n", u, v, cf_adj[u*num_nodes + v]);
                //         int curr_label = labels[v];
                //         if (curr_label < min_label)
                //         {
                //             min_v = v;
                //             min_label = curr_label;
                //         }
                //     }
                // }
                int start = edge_starts[u];
                int end = edge_starts[u + 1];
                int min_edge_idx = 0;
                for (int i = start; i < end; i++)
                {
                    int v = edge_dests[i];
                    int capacity = cf[i];
                    if (capacity > 0)
                    {
                        // printf("Residual at edge (%d, %d) is %d\n", u, v, cf_adj[u*num_nodes + v]);
                        int curr_label = labels[v];
                        if (curr_label < min_label)
                        {
                            min_v = v;
                            min_label = curr_label;
                            min_edge_idx = i;
                        }
                    }
                }
                // printf("lowest dest node is %d, with label %d, u is %d\n", min_v, min_label, u);
                if (labels[u] > min_label)
                {
                    // printf("Excess of %d is %d\n", min_v, excess[min_v]);
                    int d = min(e_prime, cf[min_edge_idx]);
                    atomicAdd(&cf[reverse_edge_index[min_edge_idx]], d);
                    atomicSub(&cf[min_edge_idx], d);
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

void Graph::globalRelabel(int num_nodes, int source, int sink, std::vector<int>& excess, std::vector<int>& labels, std::vector<int>& cf, 
                          std::vector<int>& edge_starts, std::vector<int>& edge_dests, std::vector<int>& reverse_edge_index, std::vector<bool>& marked)
{
    // ensure that source label is no more than 1 greater than dest label
    int index = 0;
    for (int u = 0; u < num_nodes; u++)
    {
        int source_label = labels[u];
        for (const Edge& e : vertices[u].outgoing_edges) {
            int v = e.dest;
            int dest_label = labels[v];
            if (source_label > dest_label + 1)
            {
                int flow = cf[index];
                excess[u] -= flow;
                excess[v] += flow;
                cf[reverse_edge_index[index]] += flow;
                cf[index] = 0;
            }
            index++;
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
        
        int start = edge_starts[u];
        int end = edge_starts[u + 1];
        for (int i = start; i < end; i++)
        {
            int v = edge_dests[i];
            if (cf[reverse_edge_index[i]] > 0 && !visited[v]) // reverse residual edge exists
            {
                labels[v] = labels[u] + 1; // 1 more than parent node height
                visited[v] = true;
                q.push(v);
            }
        }
        // for (int v = 0; v < num_nodes; v++) {
        //     if (cf_adj[v*num_nodes + u] > 0 && !visited[v]) { // reverse residual edge exists
        //         labels[v] = labels[u] + 1; // 1 more than parent node height
        //         visited[v] = true;
        //         q.push(v);
        //     }
        // }
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
    int threadsPerBlock = 256;
    int numberOfBlocks = (N / threadsPerBlock) + 1;
    // Host data structures (excess, labels, capacity/flow)
    std::vector<bool> marked(N, false);
    std::vector<int> h_excess(N, 0);
    std::vector<int> h_labels(N, 0);
    // necessary to have "flattened" ds for GPU locality
    std::vector<int> h_edge_starts(N+1, 0);
    std::vector<int> h_edge_dests(M, 0);
    std::vector<int> h_cf(M, 0); // will need to adjust this for larger test cases because we will run out of memory
    std::vector<int> h_reverse_edge_index(M, 0);
    std::unordered_map<std::pair<int,int>, int, pair_hash> edge_to_index;
    init_preflow(s);
    // printf("Done initializing preflow with excess total %d\n", excess_total);
    int index = 0;
    for (int u = 0; u < N; u++) {
        h_labels[u] = vertices[u].label;
        h_excess[u] = vertices[u].excess;
        // printf("Node %d has label %d and excess %d\n", u, h_labels[u], h_excess[u]);
        h_edge_starts[u] = index;
        for (const Edge& e : vertices[u].outgoing_edges) {
            h_edge_dests.push_back(e.dest);
            h_cf.push_back(e.capacity);
            edge_to_index[{u, e.dest}] = index;
            index++;
            // No need to add reverse flow since that should already be accounted for in init_preflow?
            // printf("Residual capacity for (%d, %d) is %d\n", e.src, e.dest, e.capacity);
        }
    }
    h_edge_starts.push_back(index); // end index of last vertex
    // For populating rev_edge_index
    for (const auto& [key, idx] : edge_to_index) {
        int u = key.first;
        int v = key.second;
        int r_idx = edge_to_index[{v, u}];
        h_reverse_edge_index[idx] = r_idx;
    }

    // Allocate device memory
    int* d_excess, *d_labels, *d_edge_starts, *d_edge_dests, *d_cf, *d_reverse_edge_index;
    cudaMalloc((void **)&d_excess, N*sizeof(int));
    cudaMalloc((void **)&d_labels, N*sizeof(int));
    cudaMalloc((void **)&d_edge_starts, (N + 1)*sizeof(int));
    cudaMalloc((void **)&d_edge_dests, M*sizeof(int));
    cudaMalloc((void **)&d_cf, M*sizeof(int));
    cudaMalloc((void **)&d_reverse_edge_index, M*sizeof(int));

    cudaMemcpy(d_excess, h_excess.data(), N*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_cf, h_cf.data(), M*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_edge_starts, h_edge_starts.data(), (N+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_edge_dests, h_edge_dests.data(), M*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_reverse_edge_index, h_reverse_edge_index.data(), M*sizeof(int), cudaMemcpyHostToDevice);

    printf("Starting loop\n");
    while (excess_total != h_excess[s] + h_excess[t])
    {
        cudaMemcpy(d_labels, h_labels.data(), N*sizeof(int), cudaMemcpyHostToDevice);
        // printf("Launching kernel\n");
        push_relabel_kernel<<<numberOfBlocks, threadsPerBlock>>>(N, s, t, d_excess, d_labels, d_cf, d_edge_starts, d_edge_dests, d_reverse_edge_index);
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
    cudaFree(d_edge_starts);
    cudaFree(d_edge_dests);
    cudaFree(d_reverse_edge_index);
    return h_excess[t];
}