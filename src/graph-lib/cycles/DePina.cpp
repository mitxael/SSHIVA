/// ******************************************************************
///	File	: DePina.cpp
///	About	: Algorithm for Minimal Cycle Basis
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// - A Faster Algorithm for Minimum Cycle Basis of Graphs (T. Kavitha1, K. Mehlhorn, D. Michail1, K. Paluch) (2006)
/// - (https://pdfs.semanticscholar.org/e4cd/f12d5ef00e0ac1fa1003e2ccbc69136ae775.pdf)
/// ******************************************************************

#include "DePina.h"													/// DePina's header file

/// Definition of variables, to be used as global elsewhere
	std::vector<std::pair<int, int>> mst;							/// minimum spanning tree (Kruskal)
	std::vector<int> ShortestPath_mst;								/// shortest path between two mst vertices
	std::vector<std::vector<ut>> S;								/// canonical basis orthogonal to linear subspace "C"
	/*std::vector<std::pair<int, double>> mst_inverse;				/// edges not belonging to minimum spanning tree
	std::vector<std::vector<int>> C;								/// linear subspace of cycles
	std::vector<std::vector<int>> S_;								/// ancillary canonical basis
	std::vector<std::vector<int>> cycleVertices;					// ordered array of vertices (to prevent duplicated cycles)*/

/// Definition of member functions

void DePina(Graph &G1) {
	
	/// Generate the Minimum Spanning Tree (mst)
	if (verbosity > 1) std::cout << "Step 1 started (generate minimum spanning tree...)" << std::endl;
	KruskalMST(G1);				
	//G1.dfsFast_caller();
	//mst = G1.DFS_spanningTree;
	//mst = { std::make_pair(0,1), std::make_pair(1,2), std::make_pair(2,4), std::make_pair(1,3), std::make_pair(3,5) };	//DePinaDE

	/*std::cout << "Spanning tree: ";
	for (int x = 0; x < mst.size(); ++x)
		std::cout << "("<< mst[x].first << "," << mst[x].second << ") ";
	std::cout << std::endl;*/

	/// Sequentially determine cycles in the MCB
	if (verbosity > 1) std::cout << "Step 2 started (generate cycles...)" << std::endl;
	int cyclesMethod = 1;	/// 0: signedGraph	1:Horton
	DePinaSelection(G1, 1);		
}

void ReferenceArrayTREE_ini(Graph &G1, bool sorted) {

	/// Initialize
	int pos_Tree = G1.E - (G1.V - 1);		/// number of edges in the tree
	int pos_NonTree = 0;					/// number of edges not in the tree

	if (!sorted)
	{
		/// Edges reference clear and resize
		G1.edgesReference.clear();				
		if(G1.edgesReference.size() != G1.E) G1.edgesReference.resize(G1.E);

		/// Create the reference array of edges, based on any Spanning Tree of G:
		/// {e1,e2,...eN € E\T} + {eN+1,eN+2,...eM € T}; both in some arbitrary but fixed order
		for (int v = 0; v < G1.V; v++) {
			for (std::list<std::pair<int, double>>::iterator it = G1.adjList[v].begin(); it != G1.adjList[v].end(); ++it) {
				/// prevent duplication of edges
				if (v <= (*it).first) {
					/// if edge doesnt belong to T put it before
					if (!edgeMST(v, (*it).first)) {
						G1.edgesReference[pos_NonTree] = std::make_tuple(v, (*it).first, (*it).second);
						pos_NonTree++;
					}
					/// if edge does belongs to T put it after
					else {
						G1.edgesReference[pos_Tree] = std::make_tuple(v, (*it).first, (*it).second);
						pos_Tree++;
					}
				}
			}
		}
	}
	else
	{
		/// create a vector containing all edges increasingly ordered by weight		
		if(G1.edgesReference.empty())
			ReferenceArray_ini(G1, true);
		std::vector<std::tuple<int, int, double>> edges = G1.edgesReference;
		std::vector<std::tuple<int, int, double>>::iterator it;

		/// Edges reference clear and resize
		G1.edgesReference.clear();				
		G1.edgesReference.resize(G1.E);
		
		/// Create the reference array of edges, based on any Spanning Tree of G:
		/// {e1,e2,...eN € E\T} + {eN+1,eN+2,...eM € T}; both in some arbitrary but fixed order
		for (it = edges.begin(); it != edges.end(); ++it)	// for all edges
		{
			/// set source as current vertex
			int src = std::get<0>(*it);
			int dst = std::get<1>(*it);
			double wgt = std::get<2>(*it);

			///prevent duplication of edges
			if (src <= dst) {
				/// if edge doesnt belong to T put it before
				if (!edgeMST(src, dst)) {
					G1.edgesReference[pos_NonTree] = std::make_tuple(src, dst, wgt);
					pos_NonTree++;
				}
				/// if edge does belongs to T put it after
				else {
					G1.edgesReference[pos_Tree] = std::make_tuple(src, dst, wgt);
					pos_Tree++;
				}
			}
		}
	}

	if (verbosity > 1)
	{
		std::cout << "Edges reference: ";
		for (std::size_t z = 0; z < G1.edgesReference.size(); ++z)
			std::cout << "(" << std::get<0>(G1.edgesReference[z]) << "," << std::get<1>(G1.edgesReference[z]) << "); ";
		std::cout << std::endl;
	}
}

void ReferenceArrayEXT_ini(Graph &G, Graph &G_extended)
{
	/// Creation of the ReferenceArray
	G_extended.edgesReference.clear();						/// Edges reference clear
	G_extended.edgesReference.resize(G_extended.E);			/// Edges reference resize
	
	/// Set limits for edges positions, based on their tree-role
	/// | 0 -> pos_tree] | pos_tree -> G1.E | G1.E -> G1.E+pos_tree | G1.E+pos_tree -> G_ext.E |
	const int V = G_extended.V;
	const int offset = V/2;
	int dim = (G_extended.E/2) - ((G_extended.V/2) - 1);
	int pos_NonTree = 0;									/// number of edges not in the tree
	int pos_Tree = pos_NonTree + dim;						/// number of edges in the tree	
	int pos_NonTree_signed = G_extended.E/2;

	/// Extend G.EdgeReference
	/// (0,4) (1,4) (1,2) (3,4)		(0,3) (0,1) (2,3) (2,4)
	/// (0,9) (1,9) (1,7) (3,9)		(0,3) (0,1) (2,3) (2,4)		(5,4) (6,4) (6,2) (8,4)		(5,8) (5,6) (7,8) (7,9)
	for(std::size_t i = 0; i < G_extended.E; ++i)
	{
		if(i < pos_Tree)
		{
			/// + to -
			std::get<0>(G_extended.edgesReference[i]) = std::get<0>(G.edgesReference[i]);
			std::get<1>(G_extended.edgesReference[i]) = std::get<1>(G.edgesReference[i]) + offset;
		}
		else
		{
			if(i < G.E)
			{
				/// + to +
				std::get<0>(G_extended.edgesReference[i]) = std::get<0>(G.edgesReference[i]);
				std::get<1>(G_extended.edgesReference[i]) = std::get<1>(G.edgesReference[i]);
			}
			else
			{
				if(i < G.E+pos_Tree)
				{
					/// - to +
					std::get<0>(G_extended.edgesReference[i]) = std::get<0>(G.edgesReference[i - G.E]) + offset;
					std::get<1>(G_extended.edgesReference[i]) = std::get<1>(G.edgesReference[i - G.E]);
				}
				else
				{
					/// - to -
					//TODO: out of bounds on G
					std::get<0>(G_extended.edgesReference[i]) = std::get<0>(G.edgesReference[i - G.E]) + offset;
					std::get<1>(G_extended.edgesReference[i]) = std::get<1>(G.edgesReference[i - G.E]) + offset;
				}
			}
		}
	}

	/// Create the reference array of edges, based on any Spanning Tree of G:
	/// {e1,e2,...eN € E\T} + {eN+1,eN+2,...eM € T}; both in some arbitrary but fixed order
	/*for (int src = 0; src < G_extended.V; src++) {
		/// For every edge in the adjacency list
		for (std::list<std::pair<int, double>>::iterator it = G_extended.adjList[src].begin(); it != G_extended.adjList[src].end(); ++it)
		{
			int dst = (*it).first;

			/// prevent duplication of edges
			if (src < dst) {
				int src_unsigned = (src < offset)? src : src-offset;
				int dst_unsigned = (dst < offset)? dst : dst-offset;
					
				/// if edge does belongs to T put it after
				if (edgeMST(src_unsigned, dst_unsigned)) {
					if (src < offset && dst < offset )
					{
						G_extended.edgesReference[pos_Tree] = std::make_pair(src, dst);
						pos_Tree++;
					}
					else {
						G_extended.edgesReference[pos_Tree_signed] = std::make_pair(src, dst);
						pos_Tree_signed++;
					}
				}
				/// if edge doesnt belong to T put it before
				else
				{
					/// G_extended.addEdge(v_pos, u_neg, w);
					if (abs(dst-src) > offset)	
					{
						G_extended.edgesReference[pos_NonTree] = std::make_pair(src, dst);
						pos_NonTree++;
					}
					/// G_extended.addEdge(v_neg, u_pos, w);
					else {
						G_extended.edgesReference[pos_NonTree_signed] = std::make_pair(src, dst);
						pos_NonTree_signed++;
					}
				}
			}
		}
	}*/
}

void witnessArray_ini(Graph &G1)
{
	/// Initialize canonical basis "S" to unit vectors, which contains vectors which are orthogonal to each other, therefore, mutually linearly independent

	/// Clear the S vector
	S.clear();					

	/// Number of cycles in minimumCycleBasis = number of edges not in MST: S[0], S[1] ... S[v]
	std::size_t N = G1.E - (G1.V - 1);	
	std::vector<ut> tmp;
	
	///for every "non-tree" edge...
	for (std::size_t i = 0; i < N; ++i)
	{
		tmp.clear();
		for (std::size_t j = 0; j < N; ++j)
		{
			if (i == j)
				tmp.push_back(1);
			else
				tmp.push_back(0);
		}
		S.push_back(tmp);
	}
}

void KruskalMST(Graph &G1)
{
	/// Algorithm based on Disjoint-sets/Union-find data structure

	/// INITIALIZE
	mst.clear();											/// clear the Minimum Spanning Tree (mst)
	int minim, maxim;
	int mst_idx = -1;										/// set MST index to a "no-yet-started" value
	std::vector<int> member(G1.V);							/// determines the belonging Set of the "i" vertex

	/// Create a vector containing all edges, sorted by weight in non-decreasing
	if(G1.edgesReference.empty())
		ReferenceArray_ini(G1, true);

	/// MAKE SET for each vertex
	for (int v = 0; v < G1.V; v++)							
		member[v] = v;

	/// DO for all edges...
	std::vector<std::tuple<int, int, double>>::iterator it;
	for (it = G1.edgesReference.begin(); it != G1.edgesReference.end() && mst_idx < G1.V; ++it)
	{
		int src = std::get<0>(*it);							/// set source as current vertex
		int dst = std::get<1>(*it);							/// set destination as the one of the current vertex
		double wgt = std::get<2>(*it);						/// set weight from source to destination

		if (member[src] != member[dst])						/// FIND if the two vertices are not members of the same Set
		{
			mst.push_back(std::make_pair(src, dst));		/// add vertex to Minimum Spanning Tree
			mst_idx++;										/// increase mst_idx to prevent useless iterations after reaching "V-1"

			if (member[src] < member[dst])					/// ensure that union is always performed into the Set with the "minimum" value
			{
				minim = member[src];
				maxim = member[dst];
			}
			else
			{
				minim = member[dst];
				maxim = member[src];
			}

			for (int j = 0; j< G1.V; ++j)					/// UNION of the two Sets (maxim e minim) into one (minim)
				if (member[j] == maxim)						/// if "j" is member of maxim, turn it into a member of minim
				member[j] = minim;
		}
	}

	if (verbosity > 1)
	{
		std::cout << "MST: ";
		for (std::size_t k = 0; k < mst.size(); k++)			//display Minimum Spanning Tree
		std::cout << "(" << mst[k].first << ", " << mst[k].second << ") ";
		std::cout << std::endl;
	}
}

bool edgeMST(int u, int v) {

	for (std::size_t i = 0; i < mst.size(); i++)
	{
		if (((mst[i].first == u) && (mst[i].second == v)) || ((mst[i].first == v) && (mst[i].second == u)))
			return true;
	}

	return false;
}

int innerProduct(std::vector<ut> &A, std::vector<ut> &B, bool bitwise, std::size_t upperbound)
{
	int tmp = 0;
	
	upperbound = (upperbound==0)? A.size() : upperbound;

	if ( (A.size() == B.size()) || (A.size() >= upperbound && A.size() >= upperbound ) )
	{
		
		for (std::size_t i = 0; i < upperbound; ++i){
			if (bitwise)
				tmp ^= A[(int)(i)] & B[(int)(i)];	// Binary : Xor (_And_)
			else
				tmp += A[(int)(i)] * B[(int)(i)];	// Decimal
		}
	}
	else
		tmp = -1;
	
	return tmp;
	//return std::inner_product(A.begin(), A.end(), B.begin(), 0);
}

void signedGraph(Graph& G1)
{
	/// Min-weight cycle can be found by shortest path computation in a two-level Graph

	/// 1. Consider an extended Graph with two copies of V: v+ [0 -> V-1] U v- [V -> 2V-1]
	const int V = G1.V * 2;
	const int E = G1.E * 2;
	const int offset = V/2;
	Graph G_extended(V);
	int dim = (E/2) - ((V/2) - 1);
	int pos_NonTree = 0;									/// number of edges not in the tree
	int pos_Tree = pos_NonTree + dim;						/// number of edges in the tree	
	int pos_NonTree_signed = E/2;
	int pos_Tree_signed = pos_NonTree_signed + dim;

	/// 2. Extend the original graph into a "signed graph":
	///		Edges that belong to S changes sides:			for e = (v, w) € S have(v0, w1) and (v1, w0)
	///		Edges that don't belong to S do not change:		for e = (v, w) !€ S have(v0, w0) and (v1, w1)
	for (std::size_t v = 0; v < offset; ++v)
	{
		std::list<std::pair<int, double>>::iterator u;
		for (u = G1.adjList[v].begin(); u != G1.adjList[v].end(); ++u)
		{
			int v_pos = static_cast<int>(v);
			int u_pos = u->first;
			double w = u->second;
			int v_neg = v_pos + offset;
			int u_neg = u_pos + offset;

			///Check vertices' index to prevent double processing
			if (v_pos < u_pos) {
				/// if edge DO NOT BELONG to Si (i.e. it's tree edge)
				if (edgeMST(v_pos, u_pos))
				{
					/// add (e1_pos, e2_pos) & (e1_neg, e2_neg) to edge set of Gi, using the same original weight
					G_extended.addEdge(v_pos, u_pos, w);
					G_extended.addEdge(v_neg, u_neg, w);
				}
				/// if edge BELONG to Si (i.e. it's cotree edge)
				else
				{
					/// add (e1_pos, e2_neg) & (e1_neg, e2_pos) to edge set of Gi, using the same original weight
					G_extended.addEdge(v_pos, u_neg, w);
					G_extended.addEdge(v_neg, u_pos, w);
				}
			}
		}
	}

	/// 3. Calculate shortest paths between v+ and v- (the change of sign ensures an ODD # of intersections with Si)
	
	/// Initialize
	spt.clear();
	spt.resize(V);
	/// Create "G1.V" shortest path trees, not "G_EXT.V"
	for(std::size_t n = 0; n < G1.V; ++n)
		dijkstra_heap(G_extended, static_cast<int>(n));
	
	/// Creation of the EdgesReference array
	ReferenceArrayEXT_ini(G1, G_extended);
	/*// Show EdgesReference array
	 *std::cout << "    ";
	for(auto it : G_extended.edgesReference)
		std::cout << "(" << it.first << "," << it.second << ") ";
	std::cout << std::endl; std::cout << "    ";
	for(auto it : G1.edgesReference)
		std::cout << "(" << it.first << "," << it.second << ") ";
	std::cout << std::endl;*/

	/// For each vertex in G
	for (std::size_t v = 0; v < G1.V; ++v)
	{
		int v_pos = static_cast<int>(v);
		int v_neg = v_pos + offset;

		/// Let pv0 be the shortest path with minimum weight among the "N" (one for each vertex v of G) paths pv (between v+ and v- in G_EXT, computed with Dijkstra)
		//#undef max	
		//double minWeight = std::numeric_limits<double>::max();

		///  Compute "V" shortest paths (v_pos->v_neg) on G_extended using Dijkstra. Total paths: V^2
		for (int root = 0; root < G1.V; ++root)
		{
			/// Determine the shortest path between "v_pos" e "v_neg"
			getShortestPath_caller(G_extended, root, v_pos, v_neg);

			//if (ShortestPathWeight_spt < minWeight)
			//{
				//minWeight = ShortestPathWeight_spt;

				if (!G_extended.cycleSpace.empty()) G_extended.removeCycleFromCS((int)G_extended.cycleSpace.size() - 1);

				//TODO: Possible improvement: avoid processing "repeated" SPT in the same loop (if spt_hash!=value)
				/// Define the new cycle
				std::vector<ut> cycle(G_extended.E);
				for (int i = 0; i < ShortestPath_spt.size(); ++i)
				{
					if (ShortestPath_spt[i] != -1)
					{
						///Conversion from ext to std
						int src = i;	//(i < offset)? i : i-offset;
						int dst = ShortestPath_spt[i];	//(ShortestPath_spt[i].first < offset)? ShortestPath_spt[i].first : ShortestPath_spt[i].first-offset;
						//G_extended.addCycleEdge(src, dst);
						cycle[G_extended.edgeReference(static_cast<int>(i), ShortestPath_spt[i])] = 1;
						//if (verbosity > 2) std::cout << ShortestPath_spt[i].first << " ";
					}
				}
				G_extended.addCycleToCS(cycle, -1);	/// -1: requires weight calculation

				/// Count occurrences of "Si" edges (with any sign)
				//TODO: Also check "odd" on edge repetitions? e.g. (u+,v+) and (u-,v-) %2 = 0 (even)
				/*int Si = 0;
				for(std::size_t idx = 0; idx < G_extended.cycleSpace.back().size(); ++idx)
				{
					int val = G_extended.cycleSpace.back()[idx];
					if( (val == 1) && (idx>pos_NonTree && idx<pos_Tree) || (idx>pos_NonTree_signed && idx<pos_Tree_signed) )
						Si++;
				}*/

				/// Determine if the cycle is: 0=even(invalid) or 1=odd(valid)
				/// e.g. 0(2,4) 4(0,3) 6(2,3) 7(0,4) Si=1 => %2 = 1
				/// (not necessary because going from v+ to v- in G_EXT changes sign so odd number is ensured)
				//if (Si%2 != 0)
				//{
					/// Move the cycle to the main graph
					std::vector<ut> newCycle(G1.E);
					for(int c = 0; c < G_extended.E; ++c)
					{
						/// Convert from G.ER to G_ext.ER
						int n = (c < G1.E)? c : c-G1.E;
						int val = (G_extended.cycleSpace[G_extended.cycleSpace.size()-1])[c];
						if(val==1)
							newCycle[n] = val;
					}
					G1.addCycleToCS(newCycle, -1); /// -1: requires weight calculation
				//}
				//else
					//G_extended.removeFromMCB(G_extended.cyclSpace.size()-1);
			//}
		}
		//std::cout << "Shortest paths computed for " << v_pos << std::endl;;
	}
	//std::cout << "All-pair shortest paths completed." << std::endl;
}

std::vector<ut> findMinCycle(Graph &G1, int e, int method)
{
	/// Compute the min-weight cycle Ci so that <Ci,Si> !=0 and Ci contains an odd number of edges from Si
	/// "Si" contains at least one edge from G\T so it forms a cycle with T.
	if (verbosity > 1) std::cout << std::endl << "    ** Analyzing edge (" << std::get<0>(G1.edgesReference[e]) << "," << std::get<1>(G1.edgesReference[e]) << "):" << std::endl;

	/// SIGNED GRAPH ALGORITHM
	if (method == 0) 
	{
		/// Generate cycles on first iteration
		//TODO: depina, depina_de
		if (e==0)
			signedGraph(G1);
	} 

	/// HORTON ALGORITHM
	else
	{
		/// Generate cycles for "e"
		//TODO: hypercube 32, 64, 128
		const bool buildER = false, buildSP = (e==0)?true:false, checkSPT = false, checkShortcuts = false, checkCommon = true, checkTiernan = false, checkIsometric = true;
		hortonGeneration(G1, G1.V, e, buildER, buildSP, checkSPT, checkShortcuts, checkCommon, checkTiernan, checkIsometric);

		/// Remove non-cycles (less than 3 edges)
		for(int i=0; i <  G1.cycleSpace.size(); ++i)
		{
			int counter = 0;
			for (std::size_t j = 0; j < G1.cycleSpace[i].size(); ++j)
				counter += G1.cycleSpace[i][j];
			if (counter < 3)
			{
				G1.removeCycleFromCS(i);
				i--;
			}
		}
	}

	/// SORT the cycles space
	G1.sortCycleSpace();
	if (verbosity > 2) G1.displayCycleSpace();

	/// FIND the first cycle (minimum weight) to be "suitable" for the MCB 
	std::size_t N = G1.E - (G1.V - 1);
	std::vector<ut>  candidate;
	for (std::size_t i = 0; i < G1.cycleSpace.size(); ++i) {
		
		/// Get the "candidate" cycle
		candidate.clear();
	
		/// Take only the first "N" edges, which are the "co-tree" to compare against Si
		std::copy(G1.cycleSpace[i].begin() + 0, G1.cycleSpace[i].begin() + N, std::back_inserter(candidate));

		/// Test if the candidate cycle is orthogonal (independent) to "S"
		/// <C[i],S[i]> != 0 : C[i] isn't null in the direction of S[i] AND has an even number of S[i] edges  (0: even; 1:odd)
		/// If the <_,_> is binary the %2 is not required (to determine the number of edges from S[j] contained in C[i])

		if ((innerProduct(candidate, S[e], true, 0) ) != 0)
		{	
			/// Put cycle in MCB
			G1.addCycleToMCB((int)i);

			if (verbosity>1)
			{
				std::cout << "    Found MinCycle for S" << e << ": ";
				for (auto& k : G1.cycleSpace[i])
					std::cout << k << ' ';
				std::cout << " *** w=" << G1.cycleWeight(static_cast<int>(i), G1.cycleSpace[i]);
				std::cout << std::endl;
				G1.displayMCB(false);
			}

			if(method==0)
				/// Remove "candidate" cycle from CycleSpace
				G1.removeCycleFromCS((int)i);		
			else
				/// Clear cycle-space for next iteration
				G1.cycleSpace.clear();	
				G1.cycleSpace_w.clear();
			
			/// found!
			break;
		}
	}

	return candidate;
}

void DePinaSelection(Graph &G1, int cyclesMethod) {

	/// INITIALIZE
	ReferenceArrayTREE_ini(G1, false);	/// reference array of edges based on MST, unsorted
	witnessArray_ini(G1);				/// witnesses							

	/// Number of cycles in minimumCycleBasis
	const std::size_t N = G1.E - (G1.V - 1);

	/// Number of cotree edges (i.e.: S[0], S[1] ... S[v])
	//const std::size_t C = G1.E - mst.size();

	/// N and cotree MUST be equal to let G have a MCB
	//assert(N > C);

	/// Do until MCB has "N" elements (i.e. it's complete)
	for (std::size_t i = 0; i < N; ++i)
	{					
		
		/// Compute the min-weight cycle Ci so that <Ci,Si>!=0 (non-orthogonal to Si, and therefore orthogonal to MCB)
		/// Si is based on "i-th" co-tree edge
		std::vector<ut>  minCycle = findMinCycle(G1, static_cast<int>(i), cyclesMethod);

		/// UPDATE_S: For all remaining (i+1 onwards) non-mst edges
		//TODO: Can be improved by relaxing the invariant to take only "k" elements from S and compute "k" cycles
		for (std::size_t j = i + 1; j < N; ++j)
		{				
			int ip = innerProduct(minCycle, S[j], true, 0);		/// 0: orthogonal; 1: non-orthogonal (iif Ci contains and odd number of edges from Sj)
			if(ip!=0)
			{
			for (std::size_t k = 0; k < N; ++k)					/// Makes Sj orthogonal to Ci, and maintains orthogonality to C1-> Ci-1
				S[j][k] ^= S[i][k];								/// Sj XOR Si : 1 iif !=
					//S[j][k] -= S[i][k] * ip;						// Sj = Sj - (Si * <Sj,Ci>) (if ip==0 then Sj remains unchanged)
																	// Saves using "if(ip != 0)" but slow down a little
			}
		}
	}

	/*if (verbosity > 2)
	{
		std::cout << "S:" << std::endl;
		for (std::size_t r = 0; r < S.size(); ++r) {
			for (std::size_t s = 0; s < S[r].size(); ++s)
				std::cout << S[r][s] << " ";
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}*/
}