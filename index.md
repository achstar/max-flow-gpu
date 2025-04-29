---
layout: home
title: "Summary"
---

Our project focuses on the implementation and parallelization of Goldberg-Tarjan's push-relabel algorithm for solving the maximum flow problem in directed flow networks. We first developed a sequential baseline implementation of push-relabel, which served as a correctness reference and a performance baseline for further optimization efforts. We then implemented and optimized a parallel version of push-relabel in CUDA, and also incorporated the global relabeling heuristic. Our optimizations of the implementation are detailed in the final report, as well as the speedup achieved for each optimization made. We are able to achieve up to 20x speedup for large graphs, and we analyze the performance benefits of each experiment we performed in our results.
