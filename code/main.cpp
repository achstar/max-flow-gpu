#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include "graph_seq.hpp"
#include "cycleTimer.h"
#include <algorithm>
#include <set>
#include <chrono>

std::string input_file;
std::string version;

void usage(const char* progname) {
    std::cerr << "Usage: " << progname << " -i <input_file> -v <seq|cuda>\n";
}

int main(int argc, char** argv)
{
    // parse command line
    int opt;
    while ((opt = getopt(argc, argv, "i:v:")) != -1) {
        switch (opt) {
        case 'i':
            input_file = optarg;
            break;
        case 'v':
            version = optarg;
            break;
        default:
            usage(argv[0]);
            return 1;
        }
    }

    if (input_file.empty() || (version != "seq" && version != "cuda" && version != "omp")) {
        usage(argv[0]);
        return 1;
    }

    std::cout << "Running " << version << " version on input: " << input_file << "\n";

    std::ifstream fin(input_file);
    if (!fin) {
        std::cerr << "Error: Cannot open input file " << input_file << std::endl;
        return 1;
    }

    std::string line;
    std::vector<Edge> edges;
    int num_nodes, num_edges, source_node, sink_node;
    while (std::getline(fin, line)) {
        if (line.empty() || line[0] == 'c') continue;

        std::istringstream iss(line);
        std::string type;
        iss >> type;

        if (type == "p") {
            std::string problem;
            iss >> problem >> num_nodes >> num_edges;
        }
        else if (type == "n") {
            int temp_node;
            std::string classifier;
            iss >> temp_node >> classifier;
            if (classifier == "s")
            {
                source_node = temp_node-1;
            }
            else if (classifier == "t")
            {
                sink_node = temp_node-1;
            }
        }
        else if (type == "a") {
            int u, v, cap;
            iss >> u >> v >> cap;
            Edge e = {.src = u-1, .dest = v-1, .capacity = cap};
            edges.push_back(e);
        }
        else
        {
            printf("Default case entered!!!\n");
        }
    }

    Graph g(num_nodes, num_edges, source_node, sink_node);
    set<pair<int, int>> edge_set;
    for (int i = 0; i < num_edges; i++)
    {
        g.addEdge(edges[i].src, edges[i].dest, edges[i].capacity, 0);
        edge_set.insert({edges[i].src, edges[i].dest});
    } 
    
    // add reverse edges
    for (int i = 0; i < num_edges; i++)
    {
        pair<int, int> reverse_edge = {edges[i].dest, edges[i].src};
        // check if reverse edge already exists
        if(edge_set.find(reverse_edge) == edge_set.end()){
            g.addEdge(edges[i].dest, edges[i].src, 0, 0);
        }
    } 
    int result;
    auto start = chrono::high_resolution_clock::now();
    double startTime = CycleTimer::currentSeconds();
    if (version == "seq")
    {
        result = g.maxFlowSeq(source_node, sink_node);
    }
    else if (version == "cuda")
    {
        result = g.maxFlowParallel(source_node, sink_node);
    }
    else if (version == "omp")
    {
        result = g.maxFlowOmp(source_node, sink_node);
    }
    double endTime = CycleTimer::currentSeconds();
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = end - start;
    printf("Time: %f\n", elapsed.count());
    printf("Result: %d\n", result);
}