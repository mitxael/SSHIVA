/// ******************************************************************
///	File	: Criccaldi.cpp
///	About	: Maximal cliques finding based on Amaldi's algorithm for MCB
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// - Efficient Deterministic Algorithms for Finding a Minimal Cycle Basis (E. Amaldi, C. Iuliano, R. Rizzi) (2010)
/// - (https://link.springer.com/chapter/10.1007/978-3-642-13036-6_30)
/// ******************************************************************

#include "Criccaldi.h"				/// Amaldi's header file

/// Definition of variables
std::vector<std::vector<int>> maxCliques;

/// Definition of methods

void Criccaldi(Graph &G1)
{
	/// Compute the Minimal Cycle Basis
	Amaldi(G1);

	/// Define cliques array
	G1.maximalCliqueSet.clear();
	std::vector<int> cliqueS_size;

	/// For each cycle in MCB
	for (int i = 0; i < G1.minimumCycleBasis.size(); ++i)
	{
		/// Count cycle edges
		std::set<int> clique;
		std::vector<int> cycleVertices(G1.V);
		int cycleE = 0;
		for (int j = 0; j < G1.minimumCycleBasis[i].size(); ++j)
		{
			if (G1.minimumCycleBasis[i][j] == 1)
			{
				cycleE++;
				cycleVertices[std::get<0>(G1.edgesReference[j])] = 1;
				cycleVertices[std::get<1>(G1.edgesReference[j])] = 1;
				clique.insert(std::get<0>(G1.edgesReference[j]));
				clique.insert(std::get<1>(G1.edgesReference[j]));
			}
		}
		int cycleV = cycleE;

		/// Density of an undirected graph
		double density = (2.0000*cycleE) / (1.0000*cycleV*(cycleV - 1.0000));
		int eMax = ((cycleV * (cycleV - 1)) / 2);
		int chordsN = eMax - cycleE;


		/// Find clique
		bool found = false;
		if (density == 1)
		{
			found = true;
		}
		else
		{
			/// Look up for chords across Edge Reference (sorted by weight in non-decreasing order)
			int foundChords = 0;
			for (std::size_t e = 0; e < G1.edgesReference.size() && foundChords < chordsN; ++e)
			{
				/// If the edge doesn't belong to the cycle
				if (G1.minimumCycleBasis[i][e] == 0)
				{
					/// Check if edge is chordal (i.e. joins two cycle's vertices)
					bool belongToCycle = (cycleVertices[std::get<0>(G1.edgesReference[e])] == 1 && cycleVertices[std::get<1>(G1.edgesReference[e])] == 1) ? true : false;

					/// if edge doesnt belong to cycle, then it's a chord
					if (belongToCycle)
					{
						clique.insert(std::get<0>(G1.edgesReference[e]));
						clique.insert(std::get<1>(G1.edgesReference[e]));
						foundChords++;
					}
				}
			}
			if (foundChords == chordsN)
				found = true;
		}
		if (found)
		{
			std::vector<int> clique_vec(clique.begin(), clique.end());
			G1.maximalCliqueSet.push_back(clique_vec);
			cliqueS_size.push_back(cycleV);
		}
	}

	/// Determine maximal among cliques
	maximalCliques(G1);

	/// Display clique
	G1.displayCliques(true, 3);
}

void maximalCliques(Graph &G1)
{
	maxCliques.clear();

	/// Look for every two cliques which union forms a bigger clique
	for (std::size_t i = 0; i < G1.maximalCliqueSet.size(); ++i)
	{
		bool expanded = false;
		for (std::size_t j = i + 1; j < G1.maximalCliqueSet.size(); ++j)
		{
			///	Determine Clique[i]-Clique[j] and Clique[j]-Clique[i]
			std::vector<int> diffSetA;
			std::set_difference(G1.maximalCliqueSet[i].begin(), G1.maximalCliqueSet[i].end(), G1.maximalCliqueSet[j].begin(), G1.maximalCliqueSet[j].end(), std::inserter(diffSetA, diffSetA.end()));
			std::vector<int> diffSetB;
			std::set_difference(G1.maximalCliqueSet[j].begin(), G1.maximalCliqueSet[j].end(), G1.maximalCliqueSet[i].begin(), G1.maximalCliqueSet[i].end(), std::inserter(diffSetB, diffSetB.end()));

			/// check if there're at least two vertices which are not connected between clique "i" and "j"
			int check = 0;
			for (int a : diffSetA)
				for (int b : diffSetB)
					if (!G1.edgeExists(a, b))
						check++;
			/// if there're no pairs of disconnected vertices between clique "i" and "j"
			if (check == 0)
			{
				/// expand clique "j" as Clique[i] U Clique[j]
				bool expanded = true;
				std::vector<int> unionAB;
				std::set_union(G1.maximalCliqueSet[i].begin(), G1.maximalCliqueSet[i].end(), G1.maximalCliqueSet[j].begin(), G1.maximalCliqueSet[j].end(), std::inserter(unionAB, unionAB.end()));
				//maxCliques.push_back(unionAB);
				G1.maximalCliqueSet[j] = unionAB;
			}
		}
		/// if there's no other clique (j) to expand the current one (i)
		if (!expanded)
			/// then the clique is -expected to be- maximal
			maxCliques.push_back(G1.maximalCliqueSet[i]);
	}

	/// Sort cliques from smallest to largest
	std::sort(maxCliques.begin(), maxCliques.end(), [](const std::vector<int> & a, const std::vector<int> & b) { return a.size() < b.size(); });

	/// Remove cliques contained in other cliques
	for (std::size_t i = 0; i < maxCliques.size(); ++i)
	{
		for (std::size_t j = maxCliques.size() - 1; j > i; --j)
		{
			std::vector<int> intersectionAB;
			std::set_intersection(maxCliques[i].begin(), maxCliques[i].end(), maxCliques[j].begin(), maxCliques[j].end(), std::inserter(intersectionAB, intersectionAB.end()));
			if (intersectionAB.size() >= maxCliques[i].size())
			{
				maxCliques.erase(maxCliques.begin() + i);
				--i;
				break;
			}
		}
	}

	G1.maximalCliqueSet = maxCliques;
}

void allCliquesBF(Graph &G1, int kmax)
{
	/// List of neighbors = G1.adjList
	auto& neighboursList = G1.adjList;
	int i = 0;

	/// at each required dimension (from low to high)
	for (int dim = 0; dim <= kmax; dim++)
	{
		if(dim==0)
		{
			for (int v1=0; v1 <= G1.V; v1++)
			{
				std::vector<int> clique = { v1 };
				G1.cliqueSpace[dim].insert(clique);
			
				///Add clique to Flitration
				int k = (int)clique.size();
				int rank = (k > 1) ? G1.lowestWeightIndex(clique) : 0;
				double weight = (k > 1) ? std::get<2>(G1.edgesReference[rank]) : 0;
				G1.PH_Filtration[rank].insert(std::make_pair(clique, weight));
			}
			
		}
		else if(dim==1)
		{
			int v1=0;
			/// for each vertex
			for (auto&& neighbours : neighboursList)
			{
				/// look for a neighbour to form an edge
				for(auto&& n : neighbours)
				{
					int v2=n.first;
					///lexicographic order
					if(v1 > v2)
					{
						std::vector<int> clique = { v1, v2 };
						G1.cliqueSpace[dim].insert(clique);

						///Add clique to Flitration
						int k = (int)clique.size();
						int rank = (k > 1) ? G1.lowestWeightIndex(clique) : 0;
						double weight = (k > 1) ? std::get<2>(G1.edgesReference[rank]) : 0;
						G1.PH_Filtration[rank].insert(std::make_pair(clique, weight));
					}
				}
				v1++;
			}
		}
		else ///(dim >= 2)
		{
			/// for each edge
			for (auto&& subclique : G1.cliqueSpace[dim-1])
			{
				std::vector<int> clique = subclique;

				/// look for a common neighbour to form a cycle 
				for (auto&& n : neighboursList[clique[0]])
				{
					int candidate = n.first;

					///consistency checks
					bool isDifferent = true;
					bool isCommon = true;
					bool isHigher = true;
					for (int k = 1; k <= dim-1; k++)
					{
						///check if already included
						if (candidate == clique[k])
						{
							isDifferent = false;
							break;
						}
						///check if it's in lexical order
						else if (candidate < clique[k])
						{
							isHigher = false;
							break;
						}
						///check if there's an edge
						else if (!G1.edgeExists(candidate, clique[k]))
						{
							isCommon = false;
							break;
						}
					}

					if (isDifferent && isCommon && isHigher)
					{
						clique.push_back(candidate);
						G1.cliqueSpace[dim].insert(clique);

						///Add clique to Flitration
						int k = dim; ///(int)clique.size();
						int rank = (k > 1) ? G1.lowestWeightIndex(clique) : 0;
						double weight = (k > 1) ? std::get<2>(G1.edgesReference[rank]) : 0;
						///ERROR HERE ?
						G1.PH_Filtration[rank].insert(std::make_pair(clique, weight));
					}
				}
			}
		}
	}
}