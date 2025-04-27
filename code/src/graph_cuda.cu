#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include "graph_seq.hpp"
#include <chrono>
#include <climits>
#include <unordered_set>
#define INF 1e9

void Graph::init_preflow_cuda(int s)
{
    vertices[s].label=N;
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

__global__ void
push_kernel(int num_nodes, int source, int sink, int* excess, int* labels, int* cf, int* edge_starts, int* edge_dests, int* reverse_edge_index, 
            int* lflag, int* v, int* edge_idx)
{
    unsigned int u = (blockIdx.x * blockDim.x) + threadIdx.x;
    if ((u < num_nodes) && (u != sink))
    {
        lflag[u] = 0;
        if ((excess[u] > 0) && (labels[u] <= num_nodes + 1))
        {
            lflag[u] = 1;
            int min_label = labels[u];
            int start = edge_starts[u];
            int end = edge_starts[u + 1];
            for (int i = start; i < end; i++)
            {
                int curr_v = edge_dests[i];
                int capacity = cf[i];
                if (capacity > 0)
                {
                    // printf("Residual at edge (%d, %d) is %d\n", u, v, cf_adj[u*num_nodes + v]);
                    int curr_label = labels[curr_v];
                    if (curr_label < min_label)
                    {
                        lflag[u] = 0;
                        v[u] = curr_v;
                        min_label = curr_label;
                        edge_idx[u] = i;
                        // printf("Index for node %d to %d is %d\n", u, curr_v, i);
                        break;
                    }
                }
            }
        }
        else
        {
            edge_idx[u] = -1;
            lflag[u] = 0;
        }
    }
}

__global__ void
color_kernel(int num_nodes, int source, int sink, int* excess, int* labels, int* cf, int* edge_starts, int* edge_dests, int* reverse_edge_index, 
               int color, int* colors, int* lflag, int* v, int* edge_idx)
{
    unsigned int u = (blockIdx.x * blockDim.x) + threadIdx.x;
    if ((u < num_nodes) && (u != sink))
    {
        int flag = lflag[u];
        int dest = v[u];
        int edge_index = edge_idx[u];
        if ((edge_index >= 0) && (colors[edge_index] == color) && (flag == 0))
        {
            // printf("Handling push for edge %d, %d of color %d, edge index %d\n", u, dest, color, edge_index);
            if (u == 6 && dest == 7)
            {
                // printf("Initial Flow values for edge 6, 7: %d, 7, 6: %d\n", cf[edge_index], cf[reverse_edge_index[edge_index]]);
            }
            int d = (excess[u] > cf[edge_index]) ? cf[edge_index] : excess[u]; // min
            cf[reverse_edge_index[edge_index]] += d;
            
            cf[edge_index] -= d;
            excess[dest] += d;
            excess[u] -= d;
            if (u == 6 && dest == 7)
            {
                // printf("Flow values for edge 6, 7: %d, 7, 6: %d\n", cf[edge_index], cf[reverse_edge_index[edge_index]]);
            }
        }
    }
}

__global__ void
relabel_kernel(int num_nodes, int sink, int* cf, int* labels, int* edge_starts, int* edge_dests, int* lflag)
{
    unsigned int u = (blockIdx.x * blockDim.x) + threadIdx.x;
    if ((u < num_nodes) && (u != sink))
    {
        int flag = lflag[u];
        if (flag == 1)
        {
            int min_label = INT_MAX;
            int start = edge_starts[u];
            int end = edge_starts[u+1];
            for (int i = start; i < end; i++)
            {
                int neighbor = edge_dests[i];
                int capacity = cf[i];
                if ((capacity > 0) && (labels[neighbor] < min_label))
                {
                    min_label = labels[neighbor];
                }
            }
            int prev_label = labels[u];
            labels[u] = min_label + 1;
            // printf("Relabeled node %d to %d from prev label %d\n", u, labels[u], prev_label);
        }
    }
}

void Graph::globalRelabel(int num_nodes, int source, int sink, std::vector<int>& excess, std::vector<int>& labels, std::vector<int>& cf, 
    std::vector<int>& edge_starts, std::vector<int>& edge_dests, std::vector<int>& reverse_edge_index, std::vector<bool>& marked)
{
    // ensure that source label is no more than 1 greater than dest label
    for (int u = 0; u < num_nodes; u++)
    {
        int source_label = labels[u];
        int start = edge_starts[u];
        int end = edge_starts[u + 1];
        for (int i = start; i < end; i++){
            int v = edge_dests[i];
            int dest_label = labels[v];
            // what does this do lol 
            if (source_label > dest_label + 1)
            {
                int flow = cf[i];
                excess[u] -= flow;
                // if (excess[u] < 0)
                // {
                //     printf("ALERT ALERT NEGATIVE NODE %d with excess %d, dest node %d\n", u, excess[u], v);
                // }
                excess[v] += flow;
                // printf("Changing excess at node %d in global relabel to %d from %d\n", v, excess[v], excess[v] - flow);
                // printf("Cf reverse residual %d with edge (%d, %d)\n", cf[reverse_edge_index[i]], u, v);
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
                if(v == source){
                    printf("global relabel hit source\n");
                }
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
        if (!visited[u]) { // if not visited and not relabeled
            // if(u == source) printf("source\n");
            // printf("Untouched node is %d, excess is %d\n", u, excess[u]);
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
    init_preflow_cuda(s);
    printf("Done initializing preflow with excess total %d\n", excess_total);
    int index = 0;
    for (int u = 0; u < N; u++) {
        h_labels[u] = vertices[u].label;
        h_excess[u] = vertices[u].excess;
        // printf("Node %d has label %d and excess %d\n", u, h_labels[u], h_excess[u]);
        h_edge_starts[u] = index;
        for (const Edge& e : vertices[u].outgoing_edges) {
            h_edge_dests[index] = e.dest;
            h_cf[index] = e.capacity;
            // printf("Edge from %d to %d has capacity %d\n", e.src, e.dest, e.capacity);
            edge_to_index[{u, e.dest}] = index;
            index++;
            // No need to add reverse flow since that should already be accounted for in init_preflow?
            // printf("Residual capacity for (%d, %d) is %d\n", e.src, e.dest, e.capacity);
        }
    }
    h_edge_starts[N] = index; // end index of last vertex
    // For populating rev_edge_index
    std::vector<int> h_colors(M, -1); // will hold color index for each edge
    std::unordered_map<int, std::unordered_set<int>> used_colors;
    std::vector<int> next_color(N, 0); // next available color for each node
    int max_color = 0;
    for (const auto& [key, idx] : edge_to_index) {
        int u = key.first;
        int v = key.second;
        int r_idx = edge_to_index[{v, u}];
        h_reverse_edge_index[idx] = r_idx;

        int color = max(next_color[u], next_color[v]);
        next_color[u] = color + 1;
        next_color[v] = color + 1;
        printf("Edge %d to %d with index %d has color %d\n", u, v, idx, color);
        h_colors[idx] = color;
        // used_colors[u].insert(color);
        // used_colors[v].insert(color);
        max_color = std::max(max_color, color);
    }
    printf("Done initializing graph colors\n");
    int num_colors = max_color + 1;
    // Allocate device memory
    int* d_excess, *d_labels, *d_edge_starts, *d_edge_dests, *d_cf, *d_reverse_edge_index, *d_colors, *d_lflag, *d_v, *d_edge_idx;
    cudaMalloc((void **)&d_excess, N*sizeof(int));
    cudaMalloc((void **)&d_labels, N*sizeof(int));
    cudaMalloc((void **)&d_edge_starts, (N + 1)*sizeof(int));
    cudaMalloc((void **)&d_edge_dests, M*sizeof(int));
    cudaMalloc((void **)&d_cf, M*sizeof(int));
    cudaMalloc((void **)&d_reverse_edge_index, M*sizeof(int));
    cudaMalloc((void **)&d_colors, M*sizeof(int));
    cudaMalloc((void **)&d_lflag, N*sizeof(int));
    cudaMalloc((void **)&d_v, N*sizeof(int));
    cudaMalloc((void **)&d_edge_idx, N*sizeof(int));

    cudaMemcpy(d_excess, h_excess.data(), N*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_cf, h_cf.data(), M*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_edge_starts, h_edge_starts.data(), (N+1)*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_edge_dests, h_edge_dests.data(), M*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_reverse_edge_index, h_reverse_edge_index.data(), M*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_colors, h_colors.data(), M*sizeof(int), cudaMemcpyHostToDevice);


    auto start = chrono::high_resolution_clock::now();
    printf("Starting loop\n");
    while (excess_total != h_excess[s] + h_excess[t])
    {
        // cudaMemcpy(d_labels, h_labels.data(), N*sizeof(int), cudaMemcpyHostToDevice);
        int cycle = N;
        while (cycle > 0)
        {
            push_kernel<<<numberOfBlocks, threadsPerBlock>>>(N, s, t, d_excess, d_labels, d_cf, d_edge_starts, d_edge_dests, d_reverse_edge_index, d_lflag, d_v, d_edge_idx);
            // cudaDeviceSynchronize(); // needed?
            for (int c = 0; c < num_colors; c++)
            {
                // printf("Relabeling color: %d\n", c);
                color_kernel<<<numberOfBlocks, threadsPerBlock>>>(N, s, t, d_excess, d_labels, d_cf, d_edge_starts, d_edge_dests, d_reverse_edge_index, c, d_colors, d_lflag, d_v, d_edge_idx);
                // cudaDeviceSynchronize(); // needed?
            }
            relabel_kernel<<<numberOfBlocks, threadsPerBlock>>>(N, t, d_cf, d_labels, d_edge_starts, d_edge_dests, d_lflag);
            // cudaDeviceSynchronize(); // needed?
            cycle--;
        }
        
        // printf("Launching kernel\n");
        // push_relabel_kernel<<<numberOfBlocks, threadsPerBlock>>>(N, s, t, d_excess, d_labels, d_cf, d_edge_starts, d_edge_dests, d_reverse_edge_index);
        // printf("Launched kernel\n");
        cudaMemcpy(h_excess.data(), d_excess, N*sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_labels.data(), d_labels, N*sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(h_cf.data(), d_cf, M*sizeof(int), cudaMemcpyDeviceToHost);
        printf("Excess total: %d\n", excess_total);
        globalRelabel(N, s, t, h_excess, h_labels, h_cf, h_edge_starts, h_edge_dests, h_reverse_edge_index, marked);
        // for(int l : h_labels){
        //     if(l >= N){
        //         // printf("setting excess_total now ....\n");
        //         excess_total = h_excess[s] + h_excess[t];
        //         break;
        //     }
        // }
        
        printf("Excess target: %d and %d\n", h_excess[t], h_excess[s]);
    }
    cudaFree(d_excess);
    cudaFree(d_labels);
    cudaFree(d_cf);
    cudaFree(d_edge_starts);
    cudaFree(d_edge_dests);
    cudaFree(d_reverse_edge_index);
    cudaFree(d_colors);
    cudaFree(d_lflag);
    cudaFree(d_v);
    cudaFree(d_edge_idx);
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    printf("Time: %f\n", elapsed.count());
    printf("Result: %d\n", h_excess[t]);
    return h_excess[t];

}