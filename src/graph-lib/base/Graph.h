/// ***********************************************************
///	File	: Graph.h
///	About	: Header of Graph class
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
/// ******************************************************************

#ifndef GRAPH_H
#define GRAPH_H

#pragma once

#include "DataStructures.h"

#include <list>						// std::list
#include <vector>					// std::vector
#include <queue>					// std:queue
#include <tuple>					// std::tuple for comparisons
#include <string>					// std::string
#include <utility>					// std::pair
#include <fstream>					// std::fstream
#include <sstream>					// std::sstream
#include <iostream>					// std::cin std::cout
#include <list>						// std::list
#include <time.h>					// time
#include <windows.h>				// file finding
#include <algorithm>				// std::max
#include <map>						// std::multimap
#include <set>						// std::set

	/// Declaration of global variables defined on ("Utils.cpp"). They will "propagate" to all files #including "Graph.h"
	extern int verbosity;

	/// Declaration of type definitions
	typedef std::uint8_t ut;									/// int from 0 to 255
	typedef unsigned short us;									/// int fom 0 to 65536
	typedef unsigned u;											/// int from 0 to 4,294,967,296
	typedef unsigned long long ull;								/// int from 0 to 18,446,744,073,709,551,616


	/// Declaration of class, and its member variables and functions
	class Graph
	{
	private: 
		Graph();												/// Default constructor
	public:
		/// ATTRIBUTES
		std::string source;										/// Source file
		int V;    												/// No. of vertices
		int E;													/// No. of edges
		int K;													/// No. of vertex of minimum clique
		double Density;											/// Density
		int Degeneracy;											/// Degeneracy
		int Cyclomatic;											/// Cyclomatic number
		int Components;											/// No. of connected components

		DataStructures::adjacencyList adjList;					/// Vector of adjacency lists (1 list per vertex)
		DataStructures::verticesTags vTags;						/// tags of vertices
		std::vector<std::tuple<int, int, double>> edgesReference;/// Reference array for edges
		std::queue<int> BFS_spanningTree;						/// Spanning tree constructed by BFS
		std::vector<std::pair<int, int>> DFS_spanningTree;		/// Spanning tree constructed by DFS
		
		std::vector<std::vector<ut>> cycleSpace;				/// Space of candidate cycles
		std::vector<double> cycleSpace_w;						/// Weights of Cycles Space
		std::vector<std::vector<ut>> minimumCycleBasis;			/// Minimum Cycle Basis
		std::vector<double> minimumCycleBasis_w;				/// Weights of Minimum Cycle Basis
		
		std::map<int, std::set<std::vector<int>>> cliqueSpace;	/// Space of all cliques
		std::vector<std::vector<int>> maximalCliqueSet;			/// Maximal Cliques Set
		std::map<std::pair<int,int>, int> weightRank;			/// Rank of weights <(u,v), rank> in non-decreasing order
		DataStructures::complexFiltration PH_Filtration; /// Filtration of simmplicial complexes
		std::map<int, DataStructures::barcodeCollection> PH_barcodes;			/// Persistent Homology data (intervals and generators)
		std::vector<std::string> algorithmS = { "Horton", "DePina", "Hybrid", "Amaldi",
				"Bron-Kerbosch (naive)", "Bron-Kerbosch (tomita)", "Bron-Kerbosch (eppstein)", "Criccaldi", ///, "Dharwadker"
				"PH (javaplex)", "PH (gudhi)", "PH (phat)", "PH (dipha)", "PH (perseus)"};
		/*std::multimap<std::string, std::string> algorithmS = {
			{"Cycles", "Horton"}, {"Cycles", "DePina"}, {"Cycles", "Hybrid"}, {"Cycles", "Amaldi"},
			{"Cliques", "Bron-Kerbosch (naive)"}, {"Cliques", "Bron-Kerbosch (tomita)"}, {"Cliques", "Bron-Kerbosch (eppstein)"}, {"Cliques", "Criccaldi"}, ///, {"Cliques", "Dharwadker"}
			{"Homology", "PH (javaplex)"}, {"Homology","PH (phat)"}, {"Homology","PH (dipha)"}, {"Homology","PH (perseus)"}, {"Homology","PH (gudhi)"} };*/

		/// METHODS
		Graph(int V);											/// Initialized constructor
		bool addEdge(int src, int dst, double w);				/// Add weighted edge to Graph
		bool removeEdge(int src, int dst);						/// Remove weighted from Graph
		bool edgeExists(int src, int dst);						/// Determine if an edge exists in G
		double edgeWeight(int src, int dst);					/// Weight of a given edge
		int edgeReference(int src, int dst);					/// Get the EdReference of a given edge
		int edgeRank(int u, int v, double w = -1);
		int lowestWeightIndex(std::vector<int> &clique);
		float density();										/// Determine the density of a Graph
		int degeneracy();										/// Degeneracy of graph
		int cyclomaticNumber();									/// Cyclomatic number of graph

		void buildGraph(boolean save2file, int n);				/// Create a Graph as an hypercube with 2^n edges
		void buildGraph(boolean save2file, int n, double density);	/// Create a graph with "n" vertices and "d" density (edges are random but uniformly distributed)
		void buildGraph(boolean save2file, int n, int m);		/// Create an euclidean graph with "n" vertices and "m" edges (edges are 100% random)
		void buildGraph(std::string file);						/// Create a Graph from a parsed file A(V,E,src,dst,w ...) or B(src,dst,w ...)
		void buildGraph();										/// Create a Graph from a file selected from menu
		void saveGraph(std::string filename, char format);		/// Export Graph to file
		void postBuild();										/// Actions taken after creating a graph
		
		void bfs(int src);										/// Breath first search
		void bfsFast();											/// Fast BFS
		bool dfs(int v, int parent, int color[], int parents[]);/// Deep first search
		void dfsFast(int v, int visited[]);						/// Fast DFS
		void dfsFast_caller();									/// Launch dfsFast
		bool isConnected(std::string mode);						/// Check if a Graph is either connected or not
		bool existCycle(void);									/// Determine if a Graph has at least one cycle
		double cycleWeight(int pos, std::vector<ut> &A);		/// Calculates the weight of a cycle
		void displayGraph();									/// Show all elements of Graph
		void displayPath(std::vector<ut> pathReference);		/// Display edges of a given path
		double displayMCB(bool verify);							/// Display all elements of the Minimum Cycle Basis
		int displayCliques(bool sorted, int Kmin);				/// Display all (maximal) cliques in the graph, osrted by size and from a given size
		void displayCycleSpace();								/// Display all elements of the Cycles Space
		double minimumCycleBasisWeight();										/// Return the value of the weight of the MCB

		void sortCycleSpace();									/// Sort the cycle space in nondecreasing order by weight
		template <typename T>
		void quickSortPlus(T * begin, T * end);					/// Order cycles [...] using an improved version of quickSort
		void quickSortCycles(int p, int r);						/// Order cycles by decreasing weight using "quick sort" algorithm
		int partition(int p, int r);							/// Ancillary quick sorting function
		void shellSortCycles(int n);							/// Order cycles by decreasing weight using "shell sort" algorithm
		void defaultSort(std::vector<std::vector<ut>>& cycle_edges, std::vector<double>& cycle_weights); /// Order cycles [...] using std::sort
		
		void addCycleToCS(std::vector<ut> &cycle, double weight); /// Add a cycle from the cycle space
		void removeCycleFromCS(int pos);						/// Remove a cycle from the cycle space
		void addHoleToMCB(std::vector<ut> &hole, double weight);/// Add a hole (1-simplex) to the minimal cycle basis
		void addCycleToMCB(int pos);							/// Add a cycle to the minimal cycle basis
		void removeCycleFromMCB(int pos);						/// Remove a cycle from the minimal cycle basis
		
		bool verifyMCB();										/// Calculate symmetric difference for all incidence vectors
		void reset();											/// Clear all ancillary arrays
		void lexicalAdjList();									/// Sort the adjacency list
		void buildWeightRank(bool nondecreasing, bool zeroBased); /// Generate a Rank of weights
		std::vector<ut> symmetricDifference_incidence(
					std::vector<ut> V1, std::vector<ut> V2);		/// Symmetric difference of two vectors<ut>
		std::vector<std::tuple<int, int, double>> symmetricDifference_edges(
					std::vector<std::tuple<int, int, double>> V1, 
					std::vector<std::tuple<int, int, double>> V2);	/// Symmetric difference of two vectors<tuple>
	};

	struct ComparePG				/// function that compares (>) the second value (weight) of two pairs (a and b)
	{
		bool operator()(const std::pair<int, double>& a, const std::pair<int, double>& b);
	};

	struct ComparePL				/// function that compares (<) the second value (weight) of two pairs (a and b)
	{
		bool operator()(const std::pair<int, double>& a, const std::pair<int, double>& b);
	};

	struct ComparePFL				/// function that compares (<) the First value of two pairs
	{
		bool operator()(const std::pair<int, double>& a, const std::pair<int, double>& b);
	};

	struct CompareTL				/// function that compares (<)the 3rd value (weight) of two tuples (a and b)
	{
		bool operator()(const std::tuple<int, int, double> &lhs, const std::tuple<int, int, double> &rhs);
	};

	struct CompareTR				/// function that compares (>)the 3rd value (weight) of two tuples (a and b)
	{
		bool operator()(const std::tuple<int, int, double> &lhs, const std::tuple<int, int, double> &rhs);
	};

	struct HashPair
	{
		size_t operator()(const std::pair<std::vector<int>, double>& p);
	};

#endif