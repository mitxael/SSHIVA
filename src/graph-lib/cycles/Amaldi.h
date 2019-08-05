/// ******************************************************************
///	File	: Amaldi.h
///	About	: Header of Amaldi class
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// - Efficient Deterministic Algorithms for Finding a Minimal Cycle Basis (E. Amaldi, C. Iuliano, R. Rizzi) (2010)
/// - (https://link.springer.com/chapter/10.1007/978-3-642-13036-6_30)
/// ******************************************************************

#ifndef AMALDI_H
#define AMALDI_H
#pragma once

#include "../base/Graph.h"				// Graph's header file
#include "Horton.h"				// Horton's header file
#include "DePina.h"				// DePina's header file
#include "Hybrid.h"				// DePina's header file
#include <vector>				// std::vector
#include <iostream>				// std::cout
#include <iterator>				// std::backinserter
#include <set>					// std::set for cliques (set of unique vertices)
//#include <list>				// std::list
//#include <tuple>				// std::tuple for Kruskal
//#include <algorithm>			// std::sort
//#include <utility>			// std::swap
//#include <numeric>			// std::inner_product


	/// Declaration of global variables defined on "DePina.cpp". They will "propagate" to all files #including "Amaldi.h"
	extern std::vector<std::pair<int, int>> mst;					/// minimum spanning tree (Kruskal)
	extern std::vector<std::pair<int, double>> ShortestPath_mst;	/// shortest path between two vertices	
	extern std::vector<std::vector<ut>> S;							/// canonical basis orthogonal to linear subspace "C"

	/// Declaration of methods
	void Amaldi(Graph &G);											/// Amaldi's method to find the minimum cycles basis
	void AdaptiveIsometricCycleExtraction(Graph &G1);				/// Sequentially select the lightest linear-independent cycles in the MCB
	void Criccaldi(Graph& G1);										/// Cliques from Amaldi

#endif