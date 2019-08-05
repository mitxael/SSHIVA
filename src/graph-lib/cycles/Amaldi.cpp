/// ******************************************************************
///	File	: Amaldi.cpp
///	About	: Amaldi's algorithm for finding the Minimum Cycle Basis
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// - Efficient Deterministic Algorithms for Finding a Minimal Cycle Basis (E. Amaldi, C. Iuliano, R. Rizzi) (2010)
/// - (https://link.springer.com/chapter/10.1007/978-3-642-13036-6_30)
/// ******************************************************************

#include "Amaldi.h"				/// Amaldi's header file

/// Definition of methods
 
void Amaldi(Graph &G1)
{
	/// [1]: Generate Horton cycles (only isometric cycles are considered) (wheel-decomposition is not checked because it's not convenient)
	if (verbosity > 1) std::cout << std::endl << "**** Step 1 started (generate cycles...)" << std::endl;
	const bool buildER = true, buildSP = true, checkSPT = true, checkShortcuts=true, checkCommon = true, checkTiernan=true, checkIsometric = true;
	hortonGeneration(G1, G1.V, G1.E, buildER, buildSP, checkSPT, checkShortcuts, checkCommon, checkTiernan, checkIsometric);

	/// [1]: Sort cycles by non-decreasing weight
	if (verbosity > 1) std::cout << std::endl << "**** Step 2 started (sort cycles...)" << std::endl;
	G1.sortCycleSpace();
	if (verbosity > 1) G1.displayCycleSpace();

	/// [2] Perform De Pina's variant for "Linear Independence Test" to iteratively build the MST and order the co-edges e1...ev (and hence the witnesses Si) adaptively
	if (verbosity > 1) std::cout << std::endl << "**** Step 3 started (select cycles...)" << std::endl;
	AdaptiveIsometricCycleExtraction(G1);
}

void AdaptiveIsometricCycleExtraction(Graph &G1)
{
	/// [2]: Number of -expected- tree edges
	const std::size_t T = G1.V - 1;

	/// [2]: Number of -expected- cotree edges = Number of cycles in minimumCycleBasis = G1.E - (G1.V - 1)
	const std::size_t N = G1.E - T;	

	/// [2]: Initialize S as empty (all zeros). S will contain the witnesses made "explicit" by an update
	S.clear();

	/// [2]: Spanning Tree edges, empty at the start
	std::vector<int> ETree(G1.E, 0);
	int ETree_elems = 0;

	/// [2]: Cotree edges, empty at the start. ECot will contain the cotress removed from S (after an update)
	std::vector<int> ECot(G1.E, 0);
	int ECot_elems = 0;

	/// ******************************************************
	/// [3]: MAKESET for each vertex for the "Union-find" construction of Spanning Tree
	int minim, maxim, cycleN = 0;
	std::vector<int> member(G1.V);
	for (int v = 0; v < G1.V; v++)	member[v] = v;
	// [3]: Do while the MCB is not complete or the CS empty
	//while( (G1.minimumCycleBasis.size() < N) && (!G1.cycleSpace.empty()) )
	/// [3] Do while all edges have not been analyzed
	while( (ETree_elems + ECot_elems < G1.E) && (!G1.cycleSpace.empty()) )
	{
		cycleN++;

		/*if (verbosity > 1)
		{
			std::cout << "     **iter#: " << cycleN << std::endl;
			std::cout << "     **Cycle: ";
			for (std::size_t c = 0; c < G1.cycleSpace[0].size(); ++c)
				std::cout << G1.cycleSpace[0][c] << " ";
			std::cout << std::endl;
			std::cout << "        Cycle details: = { ";
			for (std::size_t c = 0; c < G1.cycleSpace[0].size(); ++c) {
				if (G1.cycleSpace[0][c] == 1)
					std::cout << "(" << G1.edgesReference[c].first << "," << G1.edgesReference[c].second << ") ";
			}
			double minimumCycleBasis_single_weight = G1.cycleWeight(static_cast<int>(0), G1.cycleSpace);
			std::cout << " }  (w = " << minimumCycleBasis_single_weight << ")" << std::endl;
		}*/

		/// [4]: Extract the lightest cycle C
		const std::vector<ut>& C = G1.cycleSpace[0];
		bool cycleIsDiscarded = false;

		/// CR set of "new edges": which are in C but not in Etree nor ECot
		/// Later, CR is set to contain only "new cotree edges" (without "new tree edges")
		std::vector<int> CR(G1.E, 0);
		int CR_elems = 0;
		
		/// Add edges from C to Etree until they create a cycle	
		for (std::size_t i = 0; i < C.size(); ++i)
		{	
			/// IF cycle's edge is "new" (i.e. doesnt belong to tree nor cotree)
			if(C[i]==1 && ETree[i]==0 && ECot[i]==0)				
			{
				int src = std::get<0>(G1.edgesReference[i]);				/// set source as current vertex
				int dst = std::get<1>(G1.edgesReference[i]);				/// set destination as the one of the current vertex					
				/// FIND if the two vertices are members of the same Set
				if (member[src] == member[dst])						
				{
					/// [5]: CR = (C \ (ETree U ECot)) \ NewTreeEdges	(i.e. new cotrees)
					if(ECot[i] == 0)
					{
						CR[i] = 1;
						CR_elems++;
					}
				}
				/// if the two vertices are NOT members of the same Set
				else 												
				{
					/// [6]: ET = ET U NewTreeEdges
					ETree[i] = 1;
					ETree_elems++;
					/// ensure that union is always performed into the Set with the "minimum" value
					if (member[src] < member[dst])					
						{ minim = member[src]; maxim = member[dst]; }
					else
						{ minim = member[dst]; maxim = member[src]; }
					/// UNION of the two Sets (maxim e minim) into one (minim)
					for (int k = 0; k< G1.V; ++k)					
						if (member[k] == maxim)	member[k] = minim;	/// if "j" is member of maxim, turn it into a member of minim
					}
			}
		}
		
		/// [7]: CR(cycle's cotrees edges) is partitioned into ES and EN:

		/// [7]: ES(cotree edges belonging to some Si)
		std::vector<ut> ES(G1.E, 0);
		int ES_elems = 0;
		/// For each edge in CR
		for(std::size_t j = 0; j < CR.size(); ++j)
		{
			/// Look for that edge across all Si in S
			for(std::size_t i = 0; i < S.size(); ++i)			
				/// If edge belongs to CR and to Si
				if(CR[j] == 1 && S[i][j] == 1)
				{
					ES[j] = 1;					/// add cotree to ES
					ES_elems++;					/// increase counter of cotrees in ES
					break;						/// pass to next CR edge
				}
		}

		/// [7]: EN(cotree edges not belonging to any Si, therefore "new") : CR\ES
		std::vector<int> EN(G1.E, 0);
		/// Look in EN for not-ES edges
		for(std::size_t j = 0; j < ES.size(); ++j)
		{
			/// If edge belongs to CR but not to ES, add cotree to EN
			if(CR[j] == 1 && ES[j] == 0)	EN[j] = 1;	/// iif in CR but not in ES
			//EN[j] = CR[j] ^ ES[j];			/// XOR : 1 iif != (false case: CR0^ES1=1)
		}
		int EN_elems = CR_elems - ES_elems;		/// counter of new cotress
		
		/// DEBUG vectors CR ES EN
		/*if (verbosity == 0) {
			std::cout << "        CR_elems = " << CR_elems << " : "; for(std::size_t i = 0; i < CR.size();++i) if(CR[i]==1) std::cout << "(" << G1.edgesReference[i].first << "," << G1.edgesReference[i].second << ") "; std::cout << std::endl;
		}
		if (verbosity > 1) {
			std::cout << "        ES_elems = " << ES_elems << " : "; for(std::size_t i = 0; i < ES.size();++i) if(ES[i]==1) std::cout << "(" << G1.edgesReference[i].first << "," << G1.edgesReference[i].second << ") "; std::cout << std::endl;
		}
		if (verbosity > 1) {
			std::cout << "        EN_elems = " << EN_elems << " : "; for(std::size_t i = 0; i < EN.size();++i) if(EN[i]==1) std::cout << "(" << G1.edgesReference[i].first << "," << G1.edgesReference[i].second << ") "; std::cout << std::endl;
		}*/

		/// ***********************************************
		/// [8]: If the cycle has one or more "new cotree" edges
		if (EN_elems > 0)	
		{
			/// [9]: More than one new cotree
			if(EN_elems > 1)
			{			
				/// DEBUG info
				/*if (verbosity > 1)
					std::cout << "	The cycle has multiple cotrees in S. [Update#0]: S = S U Sj(ej,ej+1)_EN" << std::endl;*/

				/// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				/// UPDATE_#0 adding Sj = {ej, ej+1} 0 <= j < EN_elems-1
				int lastE = 0;
				for(int j = 0; j < EN_elems-1; ++j)
				{
					int n = 0;
					std::vector<ut> Sj(G1.E, 0);
					///  Put into Sj two consecutive edges of EN (j and j+1)
					for(std::size_t k = lastE; k < EN.size() && n < 2; ++k)
					{
						if(EN[k] == 1)
						{
							Sj[k] = 1;
							if(n==1) lastE = static_cast<int>(k);	/// save last index of "e" for next iter.
							n++;
						}
					}
					/// S U Sj
					S.push_back(Sj);
				}
			}

			/// [10]: Search for an "implicit" witness
			bool SjUpdated = false;		/// false : e1 is no longer used in S
			/// Look for the 1st cotree edge in EN, to use as "implicit" witness
			std::size_t e1_idx;
			for(e1_idx = 0; e1_idx < EN.size(); ++e1_idx)
			{
				/// When "e1" is found, stop and keep record of its index
				if(EN[e1_idx] == 1)
					break;
			}
			/// [10]: ES is not empty, then update explicit witnesses (using the implicit witness)
			if(ES_elems > 0)
			{
				/// DEBUG info
				/*if (verbosity > 1)
					std::cout << "	The cycle has at least one cotree in S. [Update#1]: Sj = Sj U (" << G1.edgesReference[e1_idx].first << "," << G1.edgesReference[e1_idx].second << ")."<< std::endl;*/

				/// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				/// UPDATE_#1 using an "implicit" Si (with one "1" and all "0")
				for(std::size_t j = 0; j < S.size(); ++j)
				{
					/// If ES and Sj has an odd number of common edges
					if ((innerProduct(S[j],ES, true, 0)) != 0)
					{
						/// Sj = Sj U e (i.e. add e1 to Sj)
						S[j][e1_idx] = 1;
						SjUpdated = true;	/// true : e1 is still used in S
					}
				}
			}
			
			/// [11]: Only one new cotree, therefore no witnesses needs to be updated
			/// Equivalent to selecting an "implicit" witness (with one "1" and all "0") and updating S (becomes 0)
			if(EN_elems == 1 && !SjUpdated)
			{
				/// DEBUG info
				/*if (verbosity > 1)
				 	std::cout << "	The cycle has one cotree and [Update#1] was NO performed. Cotree added: (" << G1.edgesReference[e1_idx].first << "," << G1.edgesReference[e1_idx].second << ") added to ECot." << std::endl;*/

				/// Look for the new cotree edge in EN
				/*for(int i = 0; i < EN.size(); ++i)
				{
					/// When found, add i to Cotree
					if(EN[i] == 1){
						ECot[i] = 1;
						ECot_elems++;
					}
				}*/

				/// Add e1 to ECotrees since it's no longer used in S
				ECot[e1_idx] = 1;
				ECot_elems++;
			}
		}
		/// ***********************************************
		/// [12]: If the cycle has ZERO "new cotree" edges
		else
		{
			/// Search for a presumed "explicit" Si
			int idx = -1;
			/// For all Si in S...
			for(std::size_t i = 0; i < S.size(); ++i)
			{
				/// if Si is not orthogonal to ES (i.e. have an odd number of common edges)
				if ((innerProduct(ES, S[i], true, 0)) != 0)
				{
					/// Get the index of Si
					idx = static_cast<int>(i);
					break;
				}
			}

			/// [12]: If ES is not empty (i.e. the cycle has at least one cotree into S) 
			///		  and the presumed "explicit" Si does exist (the cycle has at least one linearly independent Si)
			if(ES_elems > 0 && idx >= 0)
			{
				/// DEBUG info
				/*if (verbosity > 1) 
				{
					std::cout << "	The cycle has at least one cotree in S and it's linearly independent to S[" << idx << "]." << std::endl;
					std::cout << "	[Update#2]: S\\Si and S[j] = S[j] ^ S[" << idx << "]." << std::endl;
				}*/

				/// Get "Si" based in idx
				std::vector<ut> Si = S[idx];

				///[13]: S\Si : 1 iif is in S but not in Si
				/// UNNCESSARY, already done in [14] when j==i
				/*for(std::size_t j = 0; j < S.size(); ++j)
				{
					for(std::size_t k = 0; k < S[j].size(); ++k)
						S[j][k] ^= Si[k];			/// Sj XOR Si : 1 iif !=
				}*/

				/// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
				/// [14]: UPDATE_#2 using the "explicit" Si (with an odd number of edges on common with ES)
				for(std::size_t j = 0; j < S.size(); ++j)
				{
					/// if Si is not orthogonal to ES
					int ip = innerProduct(S[j], ES, true, 0);
					if (ip != 0)
					{
						for(std::size_t k = 0; k < S[j].size(); ++k)
							S[j][k] ^= Si[k];		/// Sj XOR Si : 1 iif !=
					}
				}

				/// [15]: For all ES edges
				/*for(std::size_t i = 0; i < ES.size(); ++i)
				{
					/// If the current ES edge is 1
					if(ES[i] == 1)
					{
						/// Put current ES edge into ECot iif it does not belong to any Sj (because removed by the update)
						if(!S.empty())
						{
							int tmp = 0;
							for(std::size_t j = 0; j < S.size(); ++j)
								//ECot[i] = ~(ECot[i] | S[j][i]);		/// NOR (1 iif 0 0)
								tmp += S[j][i];
							ECot[i] = (tmp == 0)? 1 : 0;
							if(ECot[i]==1) ECot_elems++;
						}
						else
						{
							ECot[i] = 1;
							ECot_elems++;
						}
						if (verbosity == 0 && ECot[i] == 1) std::cout << "	Cotree added: (" << G1.edgesReference[i].first << "," << G1.edgesReference[i].second << ") added to ECot." << std::endl;
					}
				}*/

				/// [15]: For all Si edges (Ensures putting into ECot all edges no longer used on S)
				for(int x=0; x < Si.size();++x)
				{
					if(Si[x] == 1)
					{
						if(!S.empty())
						{
							int tmp = 0;
							for(std::size_t j = 0; j < S.size(); ++j)
								tmp += S[j][x];
							ECot[x] = (tmp == 0)? 1 : 0;
							if(ECot[x]==1) ECot_elems++;
						}
					}
				}
			}
			/// [17]: ES is empty (the cycle has no cotree edges on S) 
			///		  or the presumed "explicit" Si DOES NOT exist (the cycle has no linearly independent Si)
			else
			{
				/// [17]: Discard cycle from cycleSpace
				G1.removeCycleFromCS(0);
				cycleIsDiscarded = true;
			}
		}

		/// [19]: At the end, if cycle is not discarded, add it into MCB
		if(!cycleIsDiscarded)
		{
			G1.addCycleToMCB(0);
			G1.removeCycleFromCS(0);

			/// DISPLAY info
			if (verbosity > 1) std::cout << "	+++++Cycle accepted." << std::endl;
		}
		else
		{
			/// DISPLAY info
			if (verbosity > 1) std::cout << "	-----Cycle discarded." << std::endl;
		}

		/// DEBUG ETree and ECotree
		#ifdef _DEBUG
		if(verbosity > 1)
		{
			std::cout << "        Tree: "; for(std::size_t i = 0; i < ETree.size();++i) if(ETree[i]==1) std::cout << "1 "; else std::cout << "0 "; std::cout << " = " << ETree_elems << std::endl;
			std::cout << "        Cotr: "; for(std::size_t i = 0; i < ECot.size();++i) if(ECot[i]==1) std::cout << "1 "; else std::cout << "0 "; std::cout << " = " << ECot_elems << std::endl;
		}
		#endif

		/// DEBUG "S" vector of witnesses
		#ifdef _DEBUG
		if(verbosity > 1)
		{
			for (std::size_t r = 0; r < S.size(); ++r) {
				std::cout << "	S[" << r <<"]: ";
				for (std::size_t s = 0; s < S[r].size(); ++s)
					std::cout << S[r][s] << " ";
				std::cout << std::endl;
			}
			std::cout << std::endl;
		}
		#endif
	}
	
	/// DEBUG summary of edges
	#ifdef _DEBUG
	if(verbosity > 1)
	{
		std::cout << "Missing cotree edges: ";
		for(int i = 0; i < ETree.size(); ++i)
			if(ETree[i] == 0 && ECot[i] == 0)
				std::cout << "(" << std::get<0>(G1.edgesReference[i]) << "," << std::get<1>(G1.edgesReference[i]) << ") ";
		std::cout << std::endl;
	}
	#endif

	/// DEBUG Spanning tree
	#ifdef _DEBUG
	if(verbosity > 1)
	{
		for(int i = 0; i < ETree.size(); ++i)
			if(ETree[i]== 1)
				std::cout <<"(" << std::get<0>(G1.edgesReference[i])<< "," << std::get<1>(G1.edgesReference[i]) << ")";
		std::cout << std::endl;
	}
	#endif
}
