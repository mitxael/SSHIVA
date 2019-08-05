/// ******************************************************************
///	File	: Horton.cpp
///	About	: Find Minimal Cycle Basis on undirected graphs
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// - A polynomial-time algorithm to find the shortest cycle basis of a graph (J. D. Horton) (1987)
/// - (https://dl.acm.org/citation.cfm?id=33351)
/// - An efficient search algorithm to find the elementary circuits of a graph (J. C. Tiernan) (1970)
/// - (https://dl.acm.org/citation.cfm?id=362819)
/// - Improved minimum cycle bases algorithms by restriction to isometric cycles (E. Amaldi, C. Iuliano, T. Jurkiewicz, K. Mehlhorn, R. Rizzi (2011)
///   (https://www.semanticscholar.org/paper/2ac65f88567af80ec7d64edc8405319cafe9e919)
/// ******************************************************************

#include "Horton.h"											/// Horton class

	/// Definition of variables, to be used as global elsewhere through "extern"
	std::vector<std::vector<std::pair<int, double>>> spt;	/// shortest path tree (found by Dijkstra)
	std::vector<int> ShortestPath_spt;						/// shortest path between two vertices
	std::vector<int> ShortestPathVertices_spt;				/// vertices involved in the shortest path between two vertices
	double ShortestPathWeight_spt;							/// value of the shortest path between two vertices
	std::tuple<int, int, int> lastSPT = {-1,-1,-1};			/// last request for determining a shortest path

/// Definition of member functions

void Horton(Graph &G1)
{
	
	/// Generate Horton's cycle space
	if (verbosity > 1) std::cout << "Step 1 started (generate cycle space...)" << std::endl;
	const bool buildER = true, buildSP = true, checkSPT = true, checkShortcuts=true, checkCommon = true, checkTiernan=true, checkIsometric = true;
	hortonGeneration(G1, G1.V, G1.E, buildER, buildSP, checkSPT, checkShortcuts, checkCommon, checkTiernan, checkIsometric);

	/// Sort cycles by non-decreasing weight
	G1.sortCycleSpace();
	if (verbosity > 1) std::cout << "Step 2 started (sort cycle space by non-decreasing weight...)" << std::endl;
	if (verbosity > 1) G1.displayCycleSpace();

	/// Select linearly-independendant cycles into the MCB
	if (verbosity > 1) std::cout << "Step 3 started (test linear independence...)" << std::endl;
	hortonSelection(G1);		
}

void ReferenceArray_ini(Graph &G1, bool sorted)
{
	
	/// initialize vectors
	G1.edgesReference.clear();
	//G1.cycleSpace.clear();
	G1.cycleSpace.resize(0);
	//G1.cycleSpace_w.clear();
	G1.cycleSpace_w.resize(0);
	
	/// create the reference array of edges
	for (int v = 0; v < G1.V; v++)
	{
		for (std::list<std::pair<int, double>>::iterator it = G1.adjList[v].begin(); it != G1.adjList[v].end(); ++it)
		{
			/// Lexicographic order by "vertex"
			if (v < (*it).first)
				G1.edgesReference.push_back(std::make_tuple(v, (*it).first, (*it).second));
			}
	}

	std::vector<std::tuple<int, int, double>> edges;
	if(sorted)
	{
		/// Lexicographic order by "weight": sort edges by weight in non-decreasing order
		std::sort(G1.edgesReference.begin(),G1.edgesReference.end(), CompareTL());
	}
}

void ShortestPaths_ini(Graph &G1)
{
	/// Reset
	spt.clear();
	spt.resize(G1.V);

	/// Generate the Shortest Path Tree (spt) from each "vertex" in G
	for (int v = 0; v < G1.V; ++v)
		dijkstra_heap(G1, v);

	/// Display all SPT
	/*for (std::size_t i = 0; i < spt.size(); ++i)
	{
		std::cout << "spt[" << i << "]: ";
		for (std::size_t j = 0; j < spt[i].size(); ++j)
			std::cout << "(" << spt[i][j].first << "," << spt[i][j].second << ") ";
		std::cout << std::endl;
	}*/
}

void dijkstra_heap(Graph &G1, int src)
{
	/// INITIALIZE
	spt[src].clear();
	#undef max															/// Prevent ambiguity between "windows.h" and "std::numeric_limits<>::max"
	const double max_double = std::numeric_limits<double>::max();
	#undef min															/// Prevent ambiguity between "windows.h" and "std::numeric_limits<>::min"
	const double min_double = std::numeric_limits<double>::epsilon();	/// OR: 0.000000000000001;
	
	/// Increase weight of each edge to reduce ties
	/*for (int a = 0; a < G1.V; ++a)
		for (std::list<std::pair<int, double>>::iterator b = G1.adjList[a].begin(); b != G1.adjList[a].end(); ++b)
			b->second += min_double;*/

	bool *visited = new bool[G1.V];								/// processing status of nodes
	int *parent = new int[G1.V];								/// parents of "visited" vertices
	double *cost = new double[G1.V];							/// min. costs from source to other nodes
	spt[src].resize(G1.V);

	for (int i = 0; i < G1.V; ++i) {							/// initialize ancillary arrays
		visited[i] = 0;											/// at the start, all nodes are unvisited
		if (i == src)											/// at the start, all costs are infinite except the "source to source"
		{
			parent[i] = -1;										/// at the start, there are no parents
			cost[i] = 0;
			spt[src][i] = std::make_pair(parent[i], cost[i]);
		}
		else
			cost[i] = max_double;								// std::numeric_limits<float>::infinity();
	}
	
	/// START ALGORITHM
	std::vector<std::pair<int, double>> heap;					/// min-heap used as a "priority queue" of all nodes, ordered increasingly by weight
	heap.push_back(std::make_pair(src, 0));						/// insert the source node (and weight=0) into the "heap" vector (as first element)
	make_heap(heap.begin(), heap.end(), ComparePG());			/// rearrange the "heap" vector to actually be a min-heap

	while (heap.size() > 0)										/// do for all nodes in the "heap", until the "heap" is empty
	{
		std::pair<int, double> p = heap[0];						/// set "p" as the first element (node, weight) of the heap
		pop_heap(heap.begin(), heap.end(), ComparePG());			/// after the insertion, rearrange the "heap" to preserve its min-heap properties
		heap.pop_back();										/// remove the last element in the heap

		if (!visited[p.first])									/// if the front node is not yet "visited"
		{
			visited[p.first] = true;							/// mark node as "visited"
			std::list<std::pair<int, double>>::iterator it;		/// define iterator for adjacent list
			for (it = G1.adjList[p.first].begin(); it != G1.adjList[p.first].end(); ++it)	/// for all nodes adjacent to "p.first", update their costs
			{
				if ((cost[p.first] + it->second) < (cost[it->first]))					/// "RELAX" if new cost is less than existing cost
				{
					//cost[it->first] = cost[p.first] + it->second;			/// update cost with the new cost.
					cost[it->first] = cost[p.first] + it->second + min_double;			/// update cost with the new cost. Add "min_double" to make "long" paths "heavier"
					parent[it->first] = p.first;										/// keep track of the parent node
					spt[src][it->first] = std::make_pair(parent[it->first], cost[it->first]);

					heap.push_back(std::make_pair(it->first, it->second + p.second));	/// insert the node (and its increased weight) into the "heap" (at the last position)
					push_heap(heap.begin(), heap.end(), ComparePG());					/// after the insertion, rearrange the "heap" to preserve its min-heap properties
				}
				
				/// Tie-breaker[...]
				if ((cost[p.first] + it->second) == (cost[it->first]))
				{
					/// No tie-break: (DEFAULT) Keep the former path.
					// Do nothing

					/// Tie-break rule #1: (DETERMINISTIC) Pick shortest-length (using a sorted-heap and an extended-weight).
					// "+ min_double" already used before

					/// Tie-break rule #2: (ARBITRARY) Pick lowest-index (lexicographic order).
					if (p.first < it->first )
					{
						cost[it->first] = cost[p.first] + it->second;						/// update cost with the new cost
						parent[it->first] = p.first;										/// keep track of the parent node
						spt[src][it->first] = std::make_pair(parent[it->first], cost[it->first]);

						heap.push_back(std::make_pair(it->first, it->second + p.second));	/// insert the node (and its increased weight) into the "heap" (at the last position)
						push_heap(heap.begin(), heap.end(), ComparePG());					/// after the insertion, rearrange the "heap" to preserve its min-heap properties
					}
				}
			}
		}
	}
}

int getShortestPathFromRoot(int root, int src, int dst)
{
	/// if the path is not a loop. (u to v)
	if (src != dst)
	{
		/// Determine "next" destination
		int nextDst = (spt[root])[dst].first;
		//double weight = (spt[root])[dst].second;

		/// if not yet at root, continue	
		if ((spt[root])[dst].first != -1)
		{	
			/// Set spt_v
			ShortestPathVertices_spt[src] = (ShortestPathVertices_spt[dst] = (ShortestPathVertices_spt[nextDst] = 1));
			
			/// Continue...
			getShortestPathFromRoot(root, src, nextDst);
		}

		/// IF there's still another node after "dst" in the path
		if (nextDst != -1)
		{
			/// set "dst" to nextHop
			ShortestPath_spt[dst] = nextDst;
		}

		return 1;
	}

	/// if the path is a loop. (u to u)
	else
	{
		return 0;
	}
}

int getShortestPathFromAny(int root, int src, int dst, bool reverse)
{
	/// Retrieve the shortest path (any vertex to any other vertex) from Dijkstra's arrays
	/// spt: [																		]
	///		  0.44	0.66	0.00	0.33	0.99	0.99	0.77	1.4		1.10	0.88	(cost)
	///		  3		0		-1		2		2		6		2		1		5		0		(parent)
	///       ^		^		^		^		^		^		^		^		^		^
	///		  0		1		2		3		4		5		6		7		8		9		(child)
	///						^__root
	///
	/// p(2,9) : 9->0  0->3  3->2	2-> -1
	/// p(4,9) : 9->0  0->3  3->2  U p(2,4,rev) = 2<-4	4<- -1
	/// p(8,9) : 9->0  0->3  3->2  U p(2,8,rev) = 5<-8	6<-5	2<-6	8<- -1

	int result = -1;
	bool dstVisited = false;

	/// if the path is not a loop. (u to v)
	if (src != dst)
	{
		/// Determine "next" destination
		int nextDst = (spt[root])[dst].first;

		/// if not yet at root, continue	
		if ((spt[root])[dst].first != -1)
		{	
			/// Set spt_v
			ShortestPathVertices_spt[src] = (ShortestPathVertices_spt[dst] = (ShortestPathVertices_spt[nextDst] = 1));

			/// Set spt_w
			//ShortestPathWeight_spt = weight;
			
			/// Continue...
			result = getShortestPathFromAny(root, src, nextDst, reverse);
		}

		/// if root is reached (dst == root)
		else
		{
			///	if the root is not the source
			if(root != src)
			{
				/// set the main source as the path's ending
				ShortestPath_spt[src] = -1;

				/// continue with the path, in "REVERSE U" mode (dst -> src)
				result = getShortestPathFromAny(root, dst, src, true);				
			}
		}

		/// IF there's still another node after "dst" in the path
		if (nextDst != -1)
		{
			/// if "standard" walking mode
			if (!reverse)
			{
				//ShortestPath_spt[dst] = (ShortestPath_spt[dst].first == -1)? std::make_pair(nextHop, weight) : std::make_pair(-1, -1.00);	

				/// if "dst" has NOT been visited, then set it to nextHop
				if(ShortestPath_spt[dst] == -1)
					ShortestPath_spt[dst] = nextDst;

				/// if "dst" has been visited, then "unvisit" it (to prevent repeated vertices)
				else
				{
					dstVisited = true;
					if(ShortestPath_spt[dst] == nextDst)
						ShortestPath_spt[dst] = -1;
					if(ShortestPath_spt[nextDst] == dst)
						ShortestPath_spt[nextDst] = -1;
				}
			}
			/// if "reverse" walking mode
			else
			{
				/// do not require checking for -1 because "reverse" entries are added before "standard" (even if walked last) 
				ShortestPath_spt[nextDst] = dst;	
				//ShortestPath_spt[nextHop] = (ShortestPath_spt[nextHop].first == -1)? std::make_pair(dst, weight) : std::make_pair(-1, -1.00);
			}
		}

		if (result >= 0)
			return result;
		else
		{
			if(!reverse && !dstVisited)
				return dst;
			else
				return -1;
		}
	}

	/// if the path is a loop. (u to u)
	else
	{
		return src;
	}
}

bool getShortestPath_caller(Graph &G1, int root, int src, int dst) {

	/// Create the SPT elements IIF the current ones belong to a different request
	if ( root != std::get<0>(lastSPT) || src != std::get<1>(lastSPT) || dst != std::get<2>(lastSPT) )
	{
		/// Set current request as "lastSPT"
		std::get<0>(lastSPT) = root; std::get<1>(lastSPT) = src; std::get<2>(lastSPT) = dst;

		/// Re-initialize variables
		ShortestPath_spt.clear();
		if(ShortestPath_spt.size()!= G1.V) ShortestPath_spt.resize(G1.V);
		//ShortestPathVertices_spt.clear();
		if(ShortestPathVertices_spt.size()!= G1.V) ShortestPathVertices_spt.resize(G1.V);
		//ShortestPathWeight_spt = 0;
		for (int i = 0; i < G1.V; ++i) {
			ShortestPath_spt[i] = -1;
			ShortestPathVertices_spt[i] = 0;				
		}

		/// Get path and weight
		if(src == root)
		{
			int result = getShortestPathFromRoot(root, src, dst);
			ShortestPathWeight_spt = (spt[root])[dst].second;
			return true;
		}
		else
		{
			int result = getShortestPathFromAny(root, src, dst, false);
			if(result == root)
				ShortestPathWeight_spt = (spt[root])[dst].second;
			else
				ShortestPathWeight_spt = (spt[root][src].second - spt[root][result].second) + (spt[root][dst].second - spt[root][result].second);
			return true;
		}
		/*if (verbosity > 2)
		{
			std::cout << "            p(" << src << "," << dst<< ") is: ";
			for (std::size_t i = 0; i < ShortestPath_spt.size(); ++i)
				if (ShortestPath_spt[i].first != -1)
					std::cout << "(" << i << "," << ShortestPath_spt[i].first << ","<< ShortestPath_spt[i].second << ") ";
			//std::cout << std::endl;
		}*/

		/// Check if a path was found
		/*for (std::size_t i = 0; i < ShortestPath_spt.size(); ++i)
			if ((ShortestPath_spt[i].first != -1) && (ShortestPath_spt[i].second != -1))
				return true;*/
	}

	/// Not created, already exists
	return false;
}

int getShortestPath_absolute(Graph &G1, int src, int dst)
{
	/// Prevent ambiguity between "windows.h" and "std::numeric_limits<>::max"
	#undef max	
	
	int SPW_idx = -1;
	double SPW = std::numeric_limits<double>::max();	// OR 999999999999.99;

	/// Construct all shortest paths and compare their weights
	for (int v = 0; v < G1.V; ++v)
	{
		getShortestPath_caller(G1, v, src, dst);
		if (ShortestPathWeight_spt < SPW) {
			SPW = ShortestPathWeight_spt;
			SPW_idx = v;
		}
	}
	return SPW_idx;
}

bool edgeSPT(int r, int u, int v) {

	for (std::size_t i = 0; i < spt.size(); i++)
		if (((static_cast<int>(i) == u) && ((spt[r])[i].first == v)) || (((static_cast<int>(i) == v) && ((spt[r])[i].first == u))))
			return true;
	return false;
}

bool checkPathIntersection(Graph &G1, int r, int u, int v) {
	
	int sx_u = nextNodeFromTo(G1, r, u);
	int sx_v = nextNodeFromTo(G1, r, v);

	if (sx_u != sx_v)
	{
		//if (verbosity > 2) std::cout << std::endl << "            Vertices " << u << " and " << v << " has only " << r << " in common" << std::endl; 
		return true;
	}
	else
	{
		//if (verbosity > 2) std::cout << std::endl << "            Vertices " << u << " and " << v << " has many vertices in common" << std::endl;
		return false;
	}

	// Calculate Path1 and Path2
	/*int commonVertices = 0;
	getShortestPath_caller(G1, r, r, u);
	std::vector<int> Path1 = ShortestPathVertices_spt;
	getShortestPath_caller(G1, r, r, v);
	std::vector<int> Path2 = ShortestPathVertices_spt;

	// check for vertices (other than root) repeated in both paths
	for (std::size_t i = 0; (i < Path1.size() && (commonVertices < 2)) ; ++i)
		if ( Path1[i] != 0 && Path2[i] != 0)
			commonVertices++;
	
	// other vertex than root is repeated in both paths	
	if (commonVertices > 1) {
		if (verbosity > 2) std::cout << std::endl <<"            Vertices " << u << " and " << v << " has many vertices in common" << std::endl;
		return false;
	}
	// if no repeated vertex was found, the intersection is satisfied
	else
	{
		if (verbosity > 2) std::cout << std::endl << "            Vertices " << u << " and " << v << " has only " << r << " in common" << std::endl;
		return true;
	}*/
}

bool pathContainSmallerVertex(Graph &G1, int src, int dst) {

	/// Generate Shortest Path Tree from src
	getShortestPath_caller(G1, src, src, dst);

	/*std::cout << "The shortest path from " << src << " to " << dst << " is: ";
	for (int i = 0; i < ShortestPath_spt.size(); ++i)
	if (ShortestPath_spt[i].first != -1)
	std::cout << "(" << i << "," << ShortestPath_spt[i].first << "," << ShortestPath_spt[i].second << ") ";*/

	/// Check all vertices in P, if any of them is bigger than src
	for (std::size_t i = 0; i < (ShortestPathVertices_spt.size()); i++)

		/// if edge exists AND one of its vertices is bigger than src AND e(src, i) isn't a co-tree edge of T(src)
		if ( ((ShortestPathVertices_spt[i] == 1) && (static_cast<int>(i) < src)) ) //&& (!(edgeSPT(src, src, static_cast<int>(i)))))
		{
			if (verbosity > 2) std::cout << "**excluded path because **" << i << ">" << src << std::endl;
			
			/// Path contains at least one vertex smaller than src
			return true;
		}
	/*std::cout << std::endl;*/

	/// Path contains only vertices bigger than src
	return false;
}

bool edgeContainsVertex(Graph &G, int e, int v)
{
	if (v == std::get<0>(G.edgesReference[e]) || v == std::get<1>(G.edgesReference[e]))
		return true;

	return false;
}

bool isNotZero(int i)
{
	return (i != 0);
}

bool edgeIsShortestPath(Graph &G, int root, int u, int v)
{
	getShortestPath_caller(G, root, u, v);

	int pathVertices = (int) std::count_if(ShortestPathVertices_spt.begin(), ShortestPathVertices_spt.end(), isNotZero);

	if (pathVertices > 2)
		return false;
	else
	{
		if (pathVertices < 2)
			return true;
		else    //(pathVertices == 2)
		{
			if ((ShortestPathVertices_spt[u] == 1 && ShortestPathVertices_spt[v] == 1))
				return true;
			else
				return false;
		}
	}
}

int nextNodeFromTo(Graph &G, int x, int v)
{
	int firstV = -1;	//ShortestPath_spt[x].first;

	/// if the ShortestPathTree already exists but for another root, recreate it
	getShortestPath_caller(G, x, x, v);

	/// search the index who has "x" in ShortestPath_spt
	for (std::size_t i = 0; i < ShortestPath_spt.size(); ++i)
	{
		if (ShortestPath_spt[i] == x)
		{
			firstV = static_cast<int>(i);
			break;
		}
	}

	return firstV;
}

bool cycleIsIsometric(Graph &G, int x, int e)
{
	/// Reduce the number of candidates cycles from H to H^
	/// Only keep isometric cycles; i.e. which have for each node "v1" an edge <v2,v3> in C so that C=C(v,v1,v2). Max.#  = n*v (v: number of cycles in the MCB = m-n+1)
	/// NOTE : Isometric cycles will be duplicated at maximum "v" times.

	/// 1. If x is an endpoint of e, say e = xv, then C = epxv = Cv, e.
	if (edgeContainsVertex(G, e, x))
		return true;

	/// 2.If x is not an endpoint of e, say e = uv, and x' = sx(u) is the first node on the best path from x to u then :
	int u = std::get<0>(G.edgesReference[e]);
	int v = std::get<1>(G.edgesReference[e]);
	int x_ = nextNodeFromTo(G, x, u);	//int x_ = nextNodeFromTo(G, v, u);

	/// (a) if x = sx'(v), then Cx', e = C,
	if (x == nextNodeFromTo(G, x_, v))
		return true;

	/// (b) if x != sx'(v) and u = sv(x') then C = Cv, xx', and
	if ((x != nextNodeFromTo(G, x_, v)) && (u == nextNodeFromTo(G, v, x_)))
		return true;

	/// (c) if x 6 = sx'(v) and u 6 = sv(x') then C is not isometric
	/// if ((x == nextNodeFromTo(G, x_, v)) && (u != nextNodeFromTo(G, v, x_)))
	return false;
}

void hortonSelection(Graph &G1) {

	/// Iteratively test the "independence" of the "lightest" candidate cycles (one-by-one), and add them to "minimumCycleBasis" if satisfactory 
	int N = G1.E - G1.V + 1;
	while ( static_cast<int>(G1.minimumCycleBasis.size()) < N && (!G1.cycleSpace.empty()) )
	{
		if (verbosity > 2) {
			std::cout << std::endl << std::endl << "*** MCB cycle has " << G1.minimumCycleBasis.size() << " cycles of " << (G1.E - G1.V + 1) << std::endl;
			G1.displayCycleSpace();
		}
		
		/// Test and add cycles until "E-V+1" cycles are added (e.g. 12 - 7 + 1 = 6)
		G1.addCycleToMCB(0);
		G1.removeCycleFromCS(0);	/// Remove "MCB" cycle from CycleSpace
			
		if (verbosity > 2) {
			std::cout << "Currently testing:" << std::endl;
				
			/// Display MCB as binary
			int lastCycle = static_cast<int>(G1.minimumCycleBasis.size()) - 1;
				for (std::size_t c = 0; c < G1.minimumCycleBasis[lastCycle].size(); ++c)
					std::cout << G1.minimumCycleBasis[lastCycle][c] << " ";
				std::cout << std::endl;
				
			/// Display MCB as reference edges
			//double minimumCycleBasis_total_weight = 0;
			std::cout << "C" << lastCycle << " = { ";
			for (std::size_t c = 0; c < G1.minimumCycleBasis[lastCycle].size(); ++c)
			{
				std::cout << "(" << std::get<0>(G1.edgesReference[c]) << "," << std::get<1>(G1.edgesReference[c]) << ") ";
			}
			double minimumCycleBasis_single_weight = G1.cycleSpace_w[lastCycle];
			//minimumCycleBasis_total_weight += minimumCycleBasis_single_weight;
			std::cout << " }  (w = " << minimumCycleBasis_single_weight << ")" << std::endl;
		}

		if (!gaussianElimination_mod2(G1.minimumCycleBasis)) {
			if (verbosity > 2) std::cout << "REJECTED " << std::endl;
				G1.removeCycleFromMCB((int)G1.minimumCycleBasis.size() - 1);
		}
		else {
			if (verbosity > 2) std::cout << "ACCEPTED " << std::endl;
		}
	}
}

void hortonGeneration(Graph &G1, int vertices, int edges, bool buildER, bool buildSP, bool checkSPT, bool checkShortcut, bool checkCommon, bool checkTiernan, bool checkIsometric)
{	
	/// Horton cycles
	/// Construct candidate cycles as: Cycle(v, x, y) = (x, y) + Path(v, x) + Path(v, y)
	/// The set has a theorical "M*N" number of elements (isometric and non-isometric)

	/// INITIALIZE
	if(buildER) ReferenceArray_ini(G1, true);					/// Set the ReferenceArray, only when Horton is executed for all edges
	if(buildSP) ShortestPaths_ini(G1);							/// Set the ShortestPaths, only when Horton is executed for all edges
	int fromE = (edges == G1.E) ? 0 : edges;					/// Set edge's lower bound according to request (one-edge or all-edges)
	int toE = (edges == G1.E) ? G1.E : edges + 1;				/// Set edge's upper bound according to request (one-edge or all-edges)

	/// For every vertex in G
	for (int v = 0; v < vertices; v++) 
	{

		/// DEBUG info
		if (verbosity > 1)
		{
			std::cout << std::endl << "    Analyzing vertex " << v << " " << std::endl;

			std::cout << "    SPT" << v << ": ";
			for (int i = 0; i < vertices; i++)
				std::cout << "(" << i << "," << (spt[v])[i].first << ")";
			std::cout << std::endl;
		}
		
		/// For every edge in G
		for (std::size_t e = fromE; static_cast<int>(e) < toE; ++e)
		{
			//int e = static_cast<int>(_e);
			int e_first = std::get<0>(G1.edgesReference[e]);
			int e_second = std::get<1>(G1.edgesReference[e]);
			double e_weight = std::get<2>(G1.edgesReference[e]);

			if (verbosity > 2) std::cout << "        Analyzing edge (" << e_first << "," << e_second <<")" << std::endl;

			/// Check that edge "e" is co-tree (i.e. it isn't contained in the SPT)
			if ( !checkSPT || (!edgeSPT(v, e_first, e_second)) )
			{

				/// Required: Consider only edges which doesn't contain "v", or does but aren't shortcuts
				bool e_contains_v = edgeContainsVertex(G1, static_cast<int>(e), v);
				if ( !checkShortcut || (!e_contains_v || (e_contains_v && !edgeIsShortestPath(G1, v, e_first, e_second))) )
				{

					/// Required: Check that the vertex "v" is the only common element between paths p[v,e1] and p[v,e2]
					if ( !checkCommon || (checkPathIntersection(G1, v, e_first, e_second)) )
					{							
					
						/// Prevent duplicates: Consider Tiernan's order "§" of vertices, and construct cycles from vertex "v" only with paths containing its following vertices in "§"
						if ( !checkTiernan || (!pathContainSmallerVertex(G1, v, e_first) && !pathContainSmallerVertex(G1, v, e_second)) )
						{

							/// Check if cycle is isometric
							if ( !checkIsometric || (cycleIsIsometric(G1, v, static_cast<int>(e))) )
							{

								/// add to the Cycle Space a new empty entry (i.e. last cycle)			
								std::vector<ut> cycle(G1.E);
								double cycle_weight = 0;

								/// Construct Horton cycles as: C[v,(e1,e2)] = p(v,e1) + p(v,e2) + (e1,e2)
								std::size_t i;

								/// Add p(v,e1)
								getShortestPath_caller(G1, v, v, e_first);
								//if (verbosity > 2) std::cout << "            Add to Cycle: ";
								cycle_weight += ShortestPathWeight_spt;
								/// add all edges of path(v,x) to last cycle
								for (i = 0; i < ShortestPath_spt.size(); i++)
								{
									if (ShortestPath_spt[i] != -1)
									{
										cycle[G1.edgeReference(static_cast<int>(i), ShortestPath_spt[i])] = 1;
										//if (verbosity > 2) std::cout << ShortestPath_spt[i].first << " ";
									}
								}

								/// Add p(v,e2)
								getShortestPath_caller(G1, v, v, e_second);
								//if (verbosity > 2) std::cout << "            Add to Cycle: ";
								cycle_weight += ShortestPathWeight_spt;
								/// add all edges of path(v,y) to last cycle
								for (i = 0; i < ShortestPath_spt.size(); i++)
								{
									if (ShortestPath_spt[i] != -1)
									{
										cycle[G1.edgeReference(static_cast<int>(i), ShortestPath_spt[i])] = 1;
										/*if (verbosity > 2) {
											std::cout << ShortestPath_spt[i].first << " ";
											if (i == ShortestPath_spt.size()-1) std::cout << std::endl;
										}*/
									}
								}

								/// Add (e1,e2)
								cycle[e] = 1;
								cycle_weight += e_weight;

								/// Add cycle weight
								G1.addCycleToCS(cycle, cycle_weight);

								if (verbosity > 1) {
									std::cout << "           *ACCEPTED edge: (" << e_first << "," << e_second << ")";
									std::cout << "        CYCLE FOR MCB: ";
									for (std::size_t c = 0; c < G1.cycleSpace[G1.cycleSpace.size() - 1].size(); ++c)
										std::cout << G1.cycleSpace[G1.cycleSpace.size() - 1][c] << " ";
									std::cout << " ... " << G1.cycleSpace_w[G1.cycleSpace_w.size() - 1] << std::endl;
								}
								//if (G1.cycleWeight(static_cast<int>(G1.cycleSpace.size()) - 1, G1.cycleSpace) == 0.58)
									//std::cout << " ISOMETRIC: " << v << "," << e_first << "," << e_second << std::endl;
							}
							else {
								if (verbosity > 2) std::cout << "            REJECTED because not isometric: (" << e_first << "," << e_second << ")" << std::endl;
							}
						}
						else {
							if (verbosity > 2) std::cout << "            REJECTED because P1 or P2 contains a vertex bigger than " << v << ": (" << e_first << "," << e_second << ")" << std::endl;
						}
					} else {
						if (verbosity > 2) std::cout << "            REJECTED because multiple common vertex: (" << e_first << "," << e_second << ")" << std::endl;
					}
				} else
				{
					if (verbosity > 2) std::cout << "            REJECTED because vertex " << v << " belongs to edge (" << e_first << "," << e_second << ")" << std::endl;
				}
			}
			else {
				if (verbosity > 2) std::cout << "            REJECTED because not co-tree (i.e. it's in SPT): (" << e_first << "," << e_second << ")" << std::endl;
			}
		}
	}
}