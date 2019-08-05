/// ******************************************************************
///	File	: Dharwadker.h
///	About	: Header of Dharwadker class
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// -  New polynomial-time algorithm for finding maximal cliques in graphs (A. Dharwadker) (2006)
/// - (http://www.dharwadker.org/clique/clique.pdf)
/// ******************************************************************

#ifndef DHARWADKER_H
#define DHARWADKER_H
#pragma once

#include "../base/Graph.h"						/// Graph's header file
#include "../cycles/Horton.h"					/// Horton's header file
 #include "../cycles/Amaldi.h"

#include <iostream> 
#include <fstream> 
#include <string> 
#include <vector> 
#include <tuple>

	/// Declaration of methods
	int Dharwadker(Graph& G);
	bool removable(std::vector<int> neighbor, std::vector<int> cover);
	int max_removable(std::vector<std::vector<int> > neighbors, std::vector<int> cover);
	std::vector<int> procedure_1(std::vector<std::vector<int> > neighbors, std::vector<int> cover);
	std::vector<int> procedure_2(std::vector<std::vector<int> > neighbors, std::vector<int> cover, int k);
	int cover_size(std::vector<int> cover);

#endif