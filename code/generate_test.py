import random
import sys
import networkx as nx

def generate_sol(filepath):
    G = nx.DiGraph()
    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('c') or line.strip() == '':
                continue
            parts = line.strip().split()
            if parts[0] == 'p':
                num_nodes = int(parts[2])
            elif parts[0] == 'n':
                if parts[2] == 's':
                    source = int(parts[1])
                elif parts[2] == 't':
                    sink = int(parts[1])
            elif parts[0] == 'a':
                u, v, cap = int(parts[1]), int(parts[2]), int(parts[3])
                G.add_edge(u, v, capacity=cap)
    flow_value, flow_dict = nx.maximum_flow(G, source, sink)
    ending = '.dimacs'
    output_path = filepath.replace(ending, '.sol')
    with open(output_path, "w") as f:
        f.write("c Flow value\n")
        f.write(f"s {flow_value}\n")

def generate_random_test(num_nodes, num_edges):
    source = 1
    sink = num_nodes

    edges = set()
    filename = f"tests/test{num_nodes}nodes{num_edges}edges.dimacs"
    with open(filename, "w") as f:
        f.write(f"c Generated {num_nodes} node test\n")
        f.write(f"p max {num_nodes} {num_edges}\n")
        f.write(f"n {source} s\n")
        f.write(f"n {sink} t\n")

        while len(edges) < num_edges:
            u = random.randint(1, num_nodes - 1)
            v = random.randint(u + 1, num_nodes)
            if u != v and (u, v) not in edges:
                seed = random.randint(1,3)
                if(seed == 1):
                    cap = random.randint(500, 1000)
                elif(seed == 2):
                    cap = random.randint(50000, 100000)
                else:
                    cap = random.randint(10000, 20000) # maybe parametrize this too?
                f.write(f"a {u} {v} {cap}\n")
                edges.add((u, v))
    
    generate_sol(filename)

def generate_grid_test(gridsize):
    filename = f"tests/grid{gridsize}x{gridsize}.dimacs"
    with open(filename, "w") as f:
        f.write(f"c {gridsize}x{gridsize} Grid\n")
        f.write(f"p max {gridsize*gridsize} {2*gridsize*(gridsize-1)}\n")
        f.write(f"n 1 s\n")
        f.write(f"n {gridsize*gridsize} t\n")
        for i in range(gridsize):
            for j in range(gridsize):
                u = i * gridsize + j + 1
                if j < gridsize - 1:
                    v = u + 1
                    f.write(f"a {u} {v} {random.randint(1, gridsize)}\n")
                if i < gridsize - 1:
                    v = u + gridsize
                    f.write(f"a {u} {v} {random.randint(1, gridsize)}\n")
    generate_sol(filename)

help_message = '''
usage: generate_test.py [-h] [-g gridsize] [-n number of nodes] [-e number of edges]

Validate the wire routing output and cost matrix

optional arguments:
  -h, --help            show this help message and exit
  -g gridsize           Grid size if specifying grid test case
  -n number of nodes    For generic random test case
  -e number of edges    For generic random test case
'''

def parse_args():
    args = sys.argv
    parsed = {}
    if '-h' in args or '--help' in args:
        print (help_message)
        sys.exit(1)
    if '-g' not in args and '-n' not in args and '-e' not in args:
        print (help_message)
        sys.exit(1)
    if '-g' in args and ('-n' in args or '-e' in args):
        print (help_message)
        sys.exit(1)
    if '-n' in args and ('-e' not in args or '-g' in args):
        print("this case")
        print (help_message)
        sys.exit(1)
    if '-g' in args:
        parsed['gridsize'] = int(args[args.index('-g') + 1])
    else:
        parsed['gridsize'] = None 
    if '-n' in args and '-e' in args:
        parsed['nodes'] = int(args[args.index('-n') + 1])
        parsed['edges'] = int(args[args.index('-e') + 1])
    return parsed

def main(args):
    if args['gridsize'] != None:
        generate_grid_test(args['gridsize'])
    else:
        generate_random_test(args['nodes'], args['edges'])

if __name__ == '__main__':
    main(parse_args())
