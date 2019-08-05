/// ******************************************************************
///	File	: Hybrid.cpp
///	About	: De Pina and Horton hybrid method for finding the Minimum Cycle Basis
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// - Implementing minimum cycle basis algorithms (K. Mehlhorn, D. Michail) (2006)
/// - (https://dl.acm.org/citation.cfm?id=1216582)
/// ******************************************************************

#include "Hybrid.h"				/// Hybrid's header file

/// Definition of methods

void Hybrid(Graph &G1) {

	/// Generate a random Spanning Tree
	if (verbosity > 1) std::cout << std::endl << "**** Step 1 started (generate minimum spanning tree...)" << std::endl;
	G1.dfsFast_caller();
	mst = G1.DFS_spanningTree;
	//KruskalMST(G1);
	
	/// Generate Horton cycles (considering only the isometric ones)
	if (verbosity > 1) std::cout << std::endl << "**** Step 2 started (generate cycles...)" << std::endl;
	const bool buildER = true, buildSP = true, checkSPT = true, checkShortcuts=true, checkCommon = true, checkTiernan=true, checkIsometric = true;
	hortonGeneration(G1, G1.V, G1.E, buildER, buildSP, checkSPT, checkShortcuts, checkCommon, checkTiernan, checkIsometric);

	/// Sort cycles by non-decreasing weight
	if (verbosity > 1) std::cout << std::endl << "**** Step 2.1 started (sort cycles...)" << std::endl;
	G1.sortCycleSpace();
	/// Display "ordered" candidate cycles
	if (verbosity > 1) G1.displayCycleSpace();

	/// Sequentially select the lightest linear-independent cycles in the MCB
	if (verbosity > 1) std::cout << std::endl << "**** Step 3 started (select cycles...)" << std::endl;
	DePinaSelection_H(G1);	
}

void convertCycleSpace(Graph &G1, bool treeBased) {

	/// Create the reference array of edges, based on any Spanning Tree of G:
	/// {e1,e2,...eN € E\T} + {eN+1,eN+2,...eM € T}; both in some arbitrary but fixed order
	std::vector<std::tuple<int, int, double>> ER_old = G1.edgesReference;
	if(treeBased)
		ReferenceArrayTREE_ini(G1, true);

	/// Convert Cycle Space from LexicoGraphic to Witness
	std::vector<std::vector<ut>> temp_CS = G1.cycleSpace;
	for (std::size_t i = 0; i < ER_old.size(); ++i)
		for (std::size_t j = 0; j <  G1.edgesReference.size(); ++j)
			if ( (std::get<0>(ER_old[i]) == std::get<0>(G1.edgesReference[j])) && (std::get<1>(ER_old[i]) == std::get<1>(G1.edgesReference[j])) )
				for (std::size_t k = 0; k < temp_CS.size(); ++k)
					G1.cycleSpace[k][j] = temp_CS[k][i];

	/*if (verbosity > 2)
	{
		std::cout << "Edges reference: ";
		for (std::size_t z = 0; z < G1.edgesReference.size(); ++z)
			std::cout << "(" << G1.edgesReference[z].first << "," << G1.edgesReference[z].second << "); ";
		std::cout << std::endl;
	}*/
}

 std::vector<ut> findMinCycle_H(Graph &G1, int i) {
	/// Compute the min-weight cycle Ci so that <Ci,Si> !=0 and Ci contains an odd number of edges from Si
	/// "Si" contains at least one edge from G\T so it forms a cycle with T.

	/// Find the first cycle (minimum weight) to be "suitable" for the MCB 
	std::size_t N = G1.E - (G1.V - 1);
	std::vector<ut>  candidate;
	for (std::size_t z = 0; z < G1.cycleSpace.size(); ++z) {
		
		/// Get the "candidate" cycle
		candidate.clear();
		std::copy(G1.cycleSpace[z].begin() + 0, G1.cycleSpace[z].begin() + N, std::back_inserter(candidate));

		/// Test if the candidate cycle is orthogonal (independent) to "S"
		/// <C[i],S[i]> != 0 : C[i] isn't null in the direction of S[i] AND has an even number of S[i] edges  (0: even; 1:odd)
		/// The <_,_> is binary, so the %2 is not required to determine the number of edges from S[j] contained in C[i]
		if ((innerProduct(candidate, S[i], true, 0) ) != 0)
		{			
			/// If the candidate pass the "test", add it into the MCB
			G1.addCycleToMCB((int)z);								/// Put cycle in MCB
			G1.removeCycleFromCS((int)z);							/// Remove "candidate" cycle from CycleSpace

			/*std::cout << "Selected: ";
			for (int m = 0; m < G1.cycleSpace[z].size(); ++m)	std::cout << G1.cycleSpace[z][m];	std::cout << " w=" << G1.cycleWeight(G1, i, G1.minimumCycleBasis);*/

			/*std::cout << std::endl << "The minimum cycle basis (MCB) is: " << std::endl;
			double minimumCycleBasis_single_weight = 0;
			double minimumCycleBasis_total_weight = 0;
			for (std::size_t r = 0; r < G1.minimumCycleBasis.size(); ++r) {
			std::cout << "C" << r << " = { ";
			for (std::size_t c = 0; c < G1.minimumCycleBasis[r].size(); ++c) {
			if (G1.minimumCycleBasis[r][c] == 1)
			std::cout << "(" << G1.edgesReference[c].first << "," << G1.edgesReference[c].second << ") ";
			}
			minimumCycleBasis_single_weight = G1.cycleWeight(G1, r, G1.minimumCycleBasis);
			minimumCycleBasis_total_weight += minimumCycleBasis_single_weight;
			std::cout << " }  (w = " << minimumCycleBasis_single_weight << ")" << std::endl;
			}
			std::cout << "The total weight of the MCB is: " << minimumCycleBasis_total_weight << std::endl;*/

			break;
		}
	}
	return candidate;
}

void DePinaSelection_H(Graph &G1) {

	/// INITIALIZE
	convertCycleSpace(G1, true);					/// reference array of edges based on MST, sorted
	witnessArray_ini(G1);							/// witnesses

	/// Number of cycles in minimumCycleBasis
	const std::size_t N = G1.E - (G1.V - 1);

	/// Number of cotree edges (i.e.: S[0], S[1] ... S[v])
	//const std::size_t C = G1.E - mst.size();

	/// N and cotree MUST be equal to let G have a MCB
	//assert(N > C);

	/// Do until MCB has "N" elements (i.e. it's complete)
	for (std::size_t i = 0; i < N; ++i)
	{	
		// Compute the min-weight cycle Ci so that <Ci,Si>!=0 (non-orthogonal to Si, and therefore orthogonal to MCB)
		// Si is based on "i-th" co-tree edge
		std::vector<ut> minCycle = findMinCycle_H(G1, static_cast<int>(i));
		
		/// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		/// Update S: For all remaining (i+1 onwards) non-mst edges
		/// update witness (i.e. make S[j] orthogonal to C[i] while keeping orthogonality to C[0] -> C[i-1])
		for (std::size_t j = i + 1; j < N; ++j)
		{				
			/// 0: orthogonal; 1: non-orthogonal (iif Ci contains and odd number of edges from Sj)
			int ip = innerProduct(minCycle, S[j], true, 0);				
			if(ip!=0)
			{
				/// Makes Sj orthogonal to Ci, and maintains orthogonality with [C1 to Ci-1]
				for (std::size_t k = 0; k < N; ++k)					
					S[j][k] ^= S[i][k];		/// Sj XOR Si : 1 iif !=
					//S[j][k] -= S[i][k] * ip;						/// Sj = Sj - (Si * <Sj,Ci>) (if ip==0 then Sj remains unchanged)
																	/// Saves using "if(ip != 0)" but slow down a little
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