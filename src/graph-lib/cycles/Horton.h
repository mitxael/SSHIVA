/// ******************************************************************
///	File	: Horton.h
///	About	: Header of Horton class
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

#ifndef HORTON_H
#define HORTON_H
#pragma once

#include "../base/Graph.h"				/// Graph's header file
#include "../base/LinearIndependence.h"	/// Linear Independence test
//#include <algorithm>			/// std::make_heap, std::pop_heap, std::push_heap, std::sort_heap, set_symmetric_difference, swap
//#include <iostream>				/// std::cin std::cout
#include <vector>				/// std::vector
//#include <limits>				/// std::numeric_limits

	/// Declaration of  methods
	void Horton(Graph &G);													/// Horton's method to find the minimum cycles basis
	//std::vector<std::tuple<int, int, double>> sortEdges(Graph &G1);			/// Sort all edges by weight, in non-decreasing order
	void ReferenceArray_ini(Graph &G, bool sorted);							/// Initialize the edges' reference array and the cycle space
	void ShortestPaths_ini(Graph &G1);										/// Determine the shortest paths for all vertices in G
	/*void cyclesWeights_ini(Graph &G1);									/// Initialize the weights of cycles	*/
	void dijkstra_heap(Graph &G, int src);									/// Dijkstra with heap algorithm
	bool getShortestPath_caller(Graph &G, int root, int src, int dst);		/// Trigger the function that determines the shortest path between two nodes
	int getShortestPathFromRoot(int root, int src, int dst);				/// Determine the shortest path between the root and any node
	int getShortestPathFromAny(int root, int src, int dst, bool reverse);	/// Determine the shortest path between two nodes, based on a given root	
	int getShortestPath_absolute(Graph &G1, int src, int dst);				/// Calculate the shortest path between two nodes, from any root
	bool edgeSPT(int r, int u, int v);										/// Check if an edge (u,v) is contained in the Shortest Path Tree generated from root (r)																			
	bool checkPathIntersection(Graph &G, int r, int u, int v);				/// Evaluate if a root node is the only intersection point between two paths
	bool pathContainSmallerVertex(Graph &G1, int src, int dst);				/// Determine if a path between "src" and "dst" has a vertex with an index lower than "src"
	bool edgeIsShortestPath(Graph &G, int root, int u, int v);				/// Determine if an edge is also the shortest path between two vertices
	bool edgeContainsVertex(Graph &G, int e, int v);						/// Determine if a vertex is part of an edge
	bool isNotZero(int i);													/// Check if a value is zero
	bool edgeIsShortestPath(Graph &G, int root, int u, int v);				/// Check if an edge is also the shortest path between its vertices
	int nextNodeFromTo(Graph &G, int x, int v);								/// Determine the 2nd-node in the path from "root" to "destination"
	bool cycleIsIsometric(Graph &G, int x, int e);							/// Determine if a cycle is isometric
	void hortonGeneration(Graph &G1, int vertices, int edges,
							bool buildER, bool buildSP, bool checkSPT,
							bool checkShortcuts, bool checkCommon, 
							bool checkTiernan, bool checkIsometric);		/// Compute Horton's cycle space
	void hortonSelection(Graph &G1);										/// Select "lightest" cycles that are linearly independant (gauss elimination)

#endif