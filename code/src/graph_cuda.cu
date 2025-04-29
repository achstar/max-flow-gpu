#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cooperative_groups.h>
#include <chrono>
#include "graph_seq.hpp"

#define INF 1e9

namespace cg = cooperative_groups;
__global__ void
push_relabel_kernel(int num_nodes, int source, int sink, int* excess, int* labels, int* cf,  int* edge_starts, int* edge_dests, int* reverse_edge_index)
{
    cg::grid_group grid = cg::this_grid();
    unsigned int u_base = (blockIdx.x * blockDim.x) + threadIdx.x;
    int cycles = num_nodes;
    while (cycles > 0)
    {
        for (unsigned int u = u_base; u < num_nodes; u += blockDim.x * gridDim.x)
        {
            if ((u < num_nodes) && (u != sink))
            {
                if ((excess[u] > 0) && labels[u] <= num_nodes + 1) 
                {
                    int e_prime = excess[u];
                    int min_label = INF;
                    int min_v = -1;
                    int start = edge_starts[u];
                    int end = edge_starts[u + 1];
                    int min_edge_idx = 0;
                    for (int i = start; i < end; i++)
                    {
                        int v = edge_dests[i];
                        int capacity = cf[i];
                        if (capacity > 0)
                        {                        
                            int curr_label = labels[v];
                            if (curr_label < min_label)
                            {
                                min_v = v;
                                min_label = curr_label;
                                min_edge_idx = i;
                            }
                        }
                    }

                    if (min_v == -1)
                    {
                        labels[u] = num_nodes;
                    }
                    else
                    {
                        if (labels[u] > min_label)
                        {
                            int d = min(e_prime, cf[min_edge_idx]);
                            atomicAdd(&cf[reverse_edge_index[min_edge_idx]], d);
                            atomicSub(&cf[min_edge_idx], d);
                            atomicAdd(&excess[min_v], d);
                            atomicSub(&excess[u], d);
                        }
                        else
                        {
                            labels[u] = min_label + 1;
                        }
                    }
                }
            }
        }
        cycles--;
        grid.sync();
    }
}

void Graph::globalRelabel(int num_nodes, int source, int sink, std::vector<int>& excess, std::vector<int>& labels, std::vector<int>& cf, 
                          std::vector<int>& edge_starts, std::vector<int>& edge_dests, std::vector<int>& reverse_edge_index)
{
    // ensure that source label is no more than 1 greater than dest label
    int index = 0;
    for (int u = 0; u < num_nodes; u++)
    {
        int source_label = labels[u];
        int start = edge_starts[u];
        int end = edge_starts[u + 1];
        for (int i = start; i < end; i++){
            int v = edge_dests[i];
            int dest_label = labels[v];
            if (source_label > dest_label + 1)
            {
                int flow = cf[i];
                excess[u] -= flow;
                excess[v] += flow;
                cf[reverse_edge_index[i]] += flow;
                cf[i] = 0;
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
    }

    for (int u = 0; u < num_nodes; u++) {
        if (!visited[u] && u != source) { // if not visited and not relabeled
            excess_total -= excess[u];
            excess[u] = 0;
        }
    }
}

int Graph::maxFlowParallel(int s, int t)
{
    // Host data structures (excess, labels, capacity/flow)
    std::vector<int> h_excess(N, 0);
    std::vector<int> h_labels(N, 0);
    // necessary to have "flattened" ds for GPU locality
    std::vector<int> h_edge_starts(N+1, 0);
    std::vector<int> h_edge_dests(M, 0);
    std::vector<int> h_cf(M, 0);
    std::vector<int> h_reverse_edge_index(M, 0);
    std::unordered_map<std::pair<int,int>, int, pair_hash> edge_to_index;
    fprintf(stderr, "before init preflow\n");
    init_preflow(s);
    fprintf(stderr, "Done initializing preflow with excess total %d\n", excess_total);
    int index = 0;
    for (int u = 0; u < N; u++) {
        h_labels[u] = vertices[u].label;
        h_excess[u] = vertices[u].excess;
        h_edge_starts[u] = index;
        for (const Edge& e : vertices[u].outgoing_edges) {
            h_edge_dests[index] = e.dest;
            h_cf[index] = e.capacity;
            edge_to_index[{u, e.dest}] = index;
            index++;
            // No need to add reverse flow since that should already be accounted for in init_preflow
        }
    }
    h_edge_starts[N] = index; // end index of last vertex
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

    // Configure the GPU
    int device;
    cudaGetDevice(&device);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, device);
    dim3 num_blocks(deviceProp.multiProcessorCount);
    dim3 block_size(1024);

    auto start = chrono::high_resolution_clock::now();
    printf("Starting loop\n");
    int iter = 0;
    while (excess_total != h_excess[s] + h_excess[t])
    {
        printf("Iteration %d: excess_total = %d, e(s) + e(t) = %d + %d = %d. Source label: %d\n", iter, excess_total, h_excess[s], h_excess[t], h_excess[s]+h_excess[t], h_labels[s]);
        cudaMemcpy(d_labels, h_labels.data(), N*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_cf, h_cf.data(), M*sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(d_excess, h_excess.data(), N*sizeof(int), cudaMemcpyHostToDevice);
        void* kernel_args[] = {&N, &s, &t, &d_excess, &d_labels, &d_cf, &d_edge_starts, &d_edge_dests, &d_reverse_edge_index};
        cudaLaunchCooperativeKernel((void*)push_relabel_kernel, num_blocks, block_size, kernel_args);
        cudaMemcpy(h_excess.data(), d_excess, N*sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_labels.data(), d_labels, N*sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_cf.data(), d_cf, M*sizeof(int), cudaMemcpyDeviceToHost);

        globalRelabel(N, s, t, h_excess, h_labels, h_cf, h_edge_starts, h_edge_dests, h_reverse_edge_index);
        iter++;
    }
    cudaFree(d_excess);
    cudaFree(d_labels);
    cudaFree(d_cf);
    cudaFree(d_edge_starts);
    cudaFree(d_edge_dests);
    cudaFree(d_reverse_edge_index);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    printf("Time: %f\n", elapsed.count());
    printf("Result: %d\n", h_excess[t]);
    return h_excess[t];

}