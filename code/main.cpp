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

    if (input_file.empty() || (version != "seq" && version != "cuda")) {
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
                source_node = temp_node;
            }
            else if (classifier == "t")
            {
                sink_node = temp_node;
            }
        }
        else if (type == "a") {
            int u, v, cap;
            iss >> u >> v >> cap;
            Edge e = {.src = u, .dest = v, .capacity = cap};
            edges.push_back(e);
        }
        else
        {
            printf("Default case entered!!!\n");
        }
    }

    Graph g(num_nodes, num_edges, source_node, sink_node);
    
    for (int i = 0; i < num_edges; i++)
    {
        g.addEdge(edges[i].src, edges[i].dest, edges[i].capacity, 0);
    } 
}