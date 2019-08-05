/// ******************************************************************
///	File	: Homology.h
///	About	: Header of Homology class
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// -  
/// ******************************************************************

#ifndef HOMOLOGY_H
#define HOMOLOGY_H
#pragma once

#include "../base/Graph.h"											/// Graph's header file
#include "../base/Utils.h"											/// Utils's header file
#include "../cliques/BronKerbosch.h"								/// Bron-Kerbosch's header file
#include "../cycles/Horton.h"										/// Horton's header file
#include "../homology/Javaplex/Javaplex.h"							/// Javaplex' header file

#include "Gudhi/Gudhi.h"											/// Gudhi
#include "PHAT/Phat.h"												/// PHAT
///#include "DiPHA/Dipha.h"											/// DiPHA
///#include "Perseus/Perseus.h"										/// Perseus


/// Declaration of class, and its member variables and functions
	class Homology
	{
	public:
		/// Declaration of methods
		static void PersistentHomology(Graph &G, int phEngine, int complexType, bool onlyInfinite, int dimension, bool saveCycles);	/// Compute homology for a given Graph
		
		static void buildComplex(Graph &G, int complexType, int dimension);

		static void buildGraphComplex(Graph &G);

		static void buildMaxCliqueComplex(Graph &G);

		static void buildMinCyclesComplex(Graph &G);

		static void buildCliquesComplex(Graph &G, bool autoAlgorithm, int dimension);

		static void buildAllCliquesComplex(Graph &G, bool autoAlgorithm);	/// Build Simplicial Complex from All Cliques
		
		static int findSubCliques_v1(Graph &G, std::map<int, std::set<std::pair<std::vector<int>, double>>> &allCliques, 
			std::vector<int> &clique, int mink, bool parallelExec);
		
		static void findSubCliques_v2(Graph &G, std::map<int, std::set<std::pair<std::vector<int>, double>>> &allCliques,
			std::map<int, std::set<std::tuple<int, double, std::vector<int>>>> &discoveringCliques, int kNow, int kMin, bool parallelExec);
		
		static void showFiltration(Graph &G);

		#undef min	/// Prevent ambiguity between "windows.h" and "std::numeric_limits<>::min"
		static double bottleneckDistance(std::string file1, std::string file2, double tolerance = std::numeric_limits<double>::min());
	};
	
#endif