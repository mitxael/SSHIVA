<p align="left">
<img width="100" height="20" src="https://travis-ci.org/mitxael/graph-suite.svg?branch=master">
</p>
<p align="left">
<img width="640" height="400" src="https://github.com/mitxael/CHoleR/graph-main.png">
</p>
<p align="left">
<img width="640" height="400" src="https://github.com/mitxael/CHoleR/graph-view.png">
</p>    
<p align="left">
<img width="800" height="400" src="https://github.com/mitxael/CHoleR/graph-benchmark.png">
</p>

    ***************************************************************************************
    * Abstract:   Collection of graph algorithms for Minimal cycles basis and Maximal cliques
    * Uses:       This software has been developed using:
    *              - MVSC 14.1 (library back-end)
    *              - C++/CLI (interface wrapper)
    *              - WPF (user front-end)
    * Author:     Michael Vasquez Otazu
    * Email:      mitxael@hotmail.it
    * Demo:       https://www.youtube.com/watch?v=NY8xB3WwQDA
    * History:  
	*             V2.1
    *             - Selection of maximum dimension to be computed in PH	
    *             V2.0 - Second release
    *             - Criccaldi alg. for Maximal cliques from Amaldi's mcb. O(m^{2} n / log n)
    *             - Bron-Kerbosch (naive) alg. for All Maximal Cliques. O(3^{n/3})
    *             - Bron-Kerbosch (tomita) alg. for All Maximal Cliques.  O(3^{n/3})
    *             - Bron-Kerbosch (eppstein) alg. for All Maximal Cliques. O(d n 3^{d/3})
    *             - Multi-threaded Benchmark for cycle/clique algorithms for a set of graphs (with .csv report)
    *             - Graph plotting with cycles/cliques visualization and selection
    *             V1.0 - First release
    *             - Support for weighted and undirected graphs
    *             - Horton alg. for Minimum Cycle Basis. O(m^3 n)
    *               (many optional improvements such as Tiernan order, Isometric cycles, etc.)
    *             - De Pina alg.for Minimum Cycle Basis. O(m^{3} + m n^{2} log n)
    *               (cycles computation using either Horton-Space or Signed-Graph)
    *             - Hybrid (Horton+DePina) alg. for Minimum Cycle Basis. O(m^{2} n + m n^{2})
    *             - Amaldi alg. for Minimum Cycle Basis. O(m^{2} n / log n)
    *             - Customizable verbosity of algorithms execution (from none to step-by-step)
    *             - Dijkstra with heap alg. O(n^{2})
    *             - DFS, BFS and many other common graph algorithms
    *             - Graph plotting
    *             - Import graph from text file (adjacency lists and single edges)
    *             - Multi-format graph export to file (adjacency lists, matrix-market, etc.)
    *             - Generate density graphs with n-vertices and uniformly-distributed random edges
    *             - Generate hypercube graphs with n^2-vertices and n*(n-1)/2 edges
    *             - Generate euclidean graphs with n-vertices and m-edges (random (x,y) positions and edges)
    ********************************* START LICENSE BLOCK *********************************
    * The MIT License (MIT)
    * Copyright (C) 2018 Michael Vasquez Otazu
    *
    * Permission is hereby granted, free of charge, to any person obtaining a copy of this 
    * software and associated documentation files (the "Software"), to deal in the Software 
    * without restriction, including without limitation the rights to use, copy, modify, merge, 
    * publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons 
    * to whom the Software is furnished to do so, subject to the following conditions:
    * 
    * The above Copyright notice and this Permission Notice shall be included in all copies 
    * or substantial portions of the Software.
    ********************************** END LICENSE BLOCK **********************************
