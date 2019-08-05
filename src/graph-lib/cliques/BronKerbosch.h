/// ******************************************************************
///	File	: BronKerbosch.h
///	About	: Header of BronKerbosch class
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// - Finding all cliques of an undirected graph (C. Bron, J. Kerbosch) (1973)
/// - (http://dl.acm.org/citation.cfm?id=362367)
/// - Worst-case time complexity for all maximal cliques (E. Tomita, A. Tanaka, H. Takahashi) (2006)
///   (http://dl.acm.org/citation.cfm?id=1217600)
/// - Listing All Maximal Cliques in Large Sparse Real-World Graphs (D. Eppstein, M Löffler, D. Strash) (2013)
///   (https://dl.acm.org/citation.cfm?id=2659459)
/// ******************************************************************


#ifndef BRONKERBOSCH_H
#define BRONKERBOSCH_H
#pragma once

#include "../base/Graph.h"				// Graph's header file
//#include "../cycles/Horton.h"				// Horton's header file
//#include "../cycles/DePina.h"				// DePina's header file
//#include "../cycles/Amaldi.h"				// Amaldi's header file
#include <vector>				// std::vector
#include <iostream>				// std::cout
#include <iterator>				// std::backinserter
#include <functional>

#include <iostream>
#include <sstream>
#include <vector>

	/// Declaration of external variables is already "inherited" from "Amaldi"
	/*extern std::vector<std::pair<int, double>> mst;				/// minimum spanning tree (Kruskal)
	extern std::vector<std::pair<int, double>> ShortestPath_mst;	/// shortest path between two vertices
	extern std::vector<std::vector<int>> S;							/// canonical basis orthogonal to linear subspace "C"
	extern std::vector<std::pair<int, double>> spt;					/// shortest path tree (found by Dijkstra)
	extern std::vector<std::pair<int, double>> ShortestPath_spt;	/// shortest path between two vertices*/

	/// Declaration of methods
	void BronKerbosch(Graph &G, int variation, int kMin=0);					/// BronKerbosch algorithm to find the Maximal Cliques (with length >= Kmin)
	void BronKerbosch_naive(Graph &G1, int Kmin,
			std::vector<int> R, std::vector<int> P, std::vector<int> X);	/// Naive version of BronKerbosch algorithm
	void BronKerbosch_tomita(Graph &G1, int Kmin,
			std::vector<int> R, std::vector<int> P, std::vector<int> X);	/// Tomita's version of BronKerbosch algorithm
	void BronKerbosch_eppstein(Graph &G1, int Kmin,
			std::vector<int> R, std::vector<int> P, std::vector<int> X);	/// Eppstein's version of BronKerbosch algorithm
	void neighbours_ini(Graph &G1);											/// Initialize array of neighbours
	int pivot(Graph& G1, std::vector<int>& P, std::vector<int>& X);			/// Vertex with lowest-degeneracy in P U X
	std::vector<int> DegeneracyOrder(Graph& G1);							/// Vector of vertices of G, sorted by degeneracy in non-decreasing order

#endif