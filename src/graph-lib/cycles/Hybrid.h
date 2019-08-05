/// ******************************************************************
///	File	: Hybrid.h
///	About	: Header of Hybrid class
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// - Implementing minimum cycle basis algorithms (K. Mehlhorn, D. Michail) (2006)
/// - (https://dl.acm.org/citation.cfm?id=1216582)
/// ******************************************************************

#ifndef HYBRID_H
#define HYBRID_H
#pragma once

#include "../base/Graph.h"		// Graph's header file
#include "Horton.h"				// Horton's header file
#include "DePina.h"				// DePina's header file
#include "Amaldi.h"				// DePina's header file
#include <vector>				// std::vector
//#include <iostream>			// std::cout
//#include <iterator>			// std::backinserter
//#include <list>				// std::list
//#include <tuple>				// std::tuple for Kruskal
//#include <algorithm>			// std::sort
//#include <utility>			// std::swap
//#include <numeric>			// std::inner_product

	/// Declaration of global variables defined on "DePina.cpp". They will "propagate" to all files #including "Amaldi.h"
	//extern std::vector<std::pair<int, int>> mst;					/// minimum spanning tree (Kruskal)
	//extern std::vector<std::pair<int, double>> ShortestPath_mst;	/// shortest path between two vertices	
	//extern std::vector<std::vector<int>> S;						/// canonical basis orthogonal to linear subspace "C"

	/// Declaration of member functions
	void Hybrid(Graph &G);											/// Hybrid method to find the minimum cycles basis
	void convertCycleSpace(Graph &G1, bool treeBased);				/// Sort incidence vectors from LexicoGraphic to Witness order
	std::vector<ut> findMinCycle_H(Graph &G1, int i);				/// Pick from "H" the lightest cycle not orthogonal to Si
	void DePinaSelection_H(Graph &G1);								/// De Pina'selection of cycles into MCB

#endif