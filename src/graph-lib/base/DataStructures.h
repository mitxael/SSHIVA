/// ***********************************************************
///	File	: DataStructures.h
///	About	: Header of Data Structures class
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
/// ******************************************************************

#ifndef DATASTRUCTURES_H
#define DATASTRUCTURES_H

#pragma once

#include <list>						/// std::list
#include <map>						/// std::multimap
#include <queue>					/// std:queue
#include <set>						/// std::set
#include <string>					/// std::string
#include <tuple>					/// std::tuple
#include <unordered_map>			/// std::unordered_map
#include <unordered_set>			/// std::unordered_set
#include <vector>					/// std::vector

	/// Declaration of class, and its member variables and functions
	class DataStructures
	{
	public:
		
		typedef std::vector<std::list<std::pair<int, double>>> adjacencyList;
		typedef std::vector<std::vector<double>> adjacencyMatrix;
		typedef std::map<int, std::string> verticesTags;

		typedef std::queue<int> intQueue;
		
		typedef std::map<int, std::set<std::pair<std::vector<int>, double>>> cliqueCollection;				/// Space of all cliques
		typedef std::map<int, std::set<std::tuple<int, double, std::vector<int>>>> foundCliqueCollection;	/// Discovered cliques
		
		typedef std::multimap<std::pair<double,double>, std::vector<std::vector<int>>> barcodeCollection;	/// Barcodes of Persistent Homology
		
		/// PH_Filtration: HashTable structure to store ordered cliques as <rank, <clique, weight>>
		/*std::map<int, std::set<std::pair<std::vector<int>, double>>> allCliques;	// Unique key with multiple values (all !=)
		std::map<int, std::unordered_set<std::pair<std::vector<int>, double>>> allCliques;	// Unique key with multiple values (all != and unordered)
		std::map<std::vector<int>, std::pair<int, double>> allCliques;				// Unique key
		std::set< std::tuple<int, std::vector<int>, double> > allCliques;				// Unique tuple-key
		std::multimap<int, std::pair<std::vector<int>, double>> allCliques;			// Multiple key-value (even =)
		// Add std::less<> to map to change order*/
		typedef std::map<int, std::set<std::pair<std::vector<int>, double>>> complexFiltration;

		std::vector<std::pair<std::string, std::string>> vectorOfStringPairs;

	};

#endif
