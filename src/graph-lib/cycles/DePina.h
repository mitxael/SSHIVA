/// ******************************************************************
///	File	: DePina.h
///	About	: Header of DePina class
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// - A Faster Algorithm for Minimum Cycle Basis of Graphs (T. Kavitha1, K. Mehlhorn, D. Michail1, K. Paluch) (2006)
/// - (https://pdfs.semanticscholar.org/e4cd/f12d5ef00e0ac1fa1003e2ccbc69136ae775.pdf)
/// ******************************************************************

#ifndef DEPINA_H
#define DEPINA_H
#pragma once

	#include "../base/Graph.h"				// Graph's header file
	#include "Horton.h"				// Horton's header file
	#include <vector>				// std::vector
	//#include <list>					// std::list
	//#include <tuple>				// std::tuple for Kruskal
	//#include <algorithm>			// std::sort
	//#include <iostream>				// std::cout
	//#include <numeric>				// std::inner_product
	//#include <iterator>				// std::backinserter

	/// Declaration of global variables defined on "Horton.cpp". They will "propagate" to all files #including "DePina.h"
	extern std::vector<std::vector<std::pair<int, double>>> spt;				/// shortest path tree (found by Dijkstra)
	extern std::vector<int> ShortestPath_spt;									/// shortest path between two vertices
	extern std::vector<int> ShortestPathVertices_spt;							/// Vertices involved in the shortest path between two vertices
	extern double ShortestPathWeight_spt;										/// Value of the shortest path between two vertices

	/// Declaration of methods
	void DePina(Graph &G);														/// De Pina's method to find the minimum cycles basis
	void KruskalMST(Graph &G);													/// Kruskal algorithm to generate a Minimum Spanning Tree
	int innerProduct(std::vector<ut> &A, std::vector<ut> &B, bool bitwise, std::size_t upbound);		/// Inner product of two vectors with the same dimension
	bool edgeMST(int u, int v);													/// Check if an edge (u,v) is contained in the Minimum Spanning Tree
	void ReferenceArrayTREE_ini(Graph &G1, bool sorted);						/// Create the reference array of edges, based on a spanning tree
	void ReferenceArrayEXT_ini(Graph &G, Graph &G_extended);					/// Create the extended reference array of edges
	void witnessArray_ini(Graph &G1);											/// Initilize the array of witness (Si)
	void signedGraph(Graph& G1);												/// Extended graph
	std::vector<ut> findMinCycle(Graph &G1, int i, int method);					/// Compute the min-weight cycle Ci, based on the "i-th" cotree edge
	void DePinaSelection(Graph &G1, int cyclesMethod);							/// Sequentially select minimum linearly-independent cycles into MCB

#endif