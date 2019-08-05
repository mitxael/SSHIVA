/// ******************************************************************
///	File	: BronKerbosch.cpp
///	About	: BronKerbosch's method for finding the Maximal Cliques
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// - Finding all cliques of an undirected graph (C. Bron, J. Kerbosch) (1973)
/// - (http://dl.acm.org/citation.cfm?id=362367)
/// - Worst-case time complexity for all maximal cliques (E. Tomita, A. Tanaka, H. Takahashi) (2006)
///   (http://dl.acm.org/citation.cfm?id=1217600)
/// - Listing All Maximal Cliques in Large Sparse Real-World Graphs (D. Eppstein, M Löffler, D. Strash) (2013)
///   (https://dl.acm.org/citation.cfm?id=2659459)
/// ******************************************************************

#include "BronKerbosch.h"													/// BronKerbosch's header file

	/// Declaration of member variables
	std::vector<std::vector<int>> neighbours;								/// []adjacent_vertices of each vertex in G

void BronKerbosch(Graph& G1, int variation, int kMin)
{	
	/// Prepare parameters
	if(kMin == 0) kMin = G1.K;
	else G1.K = kMin;
	neighbours_ini(G1);
	std::vector<int> R, P, X;
	for (int i = 0; i < G1.V; i++)
		P.push_back(i);
	
		/// Execute algorithm 
		switch (variation)
		{
			case 0:
				BronKerbosch_naive(G1, kMin, R, P, X);
				break;
			case 1:
				BronKerbosch_tomita(G1, kMin, R, P, X);
				break;
			case 2:
				BronKerbosch_eppstein(G1, kMin, R, P, X);
				break;
			default:
				std::cout << "Invalid option. Bron-Kerbosch has only 3 variants." << std::endl;
				break;
		}
	/// Display cliqueS
	bool sorted = true;
	int k = kMin;
	G1.displayCliques(sorted, kMin);
}

void neighbours_ini(Graph &G1)
{
	/// Create adjacency list without weights
	neighbours.clear();
	neighbours.resize(G1.V);
	for (int v = 0; v < G1.V; v++)
	{
		for (std::list<std::pair<int, double>>::iterator it = G1.adjList[v].begin(); it != G1.adjList[v].end(); ++it)
		{
			neighbours[v].push_back((*it).first);
		}
	}
}

void BronKerbosch_naive(Graph &G1, int Kmin, std::vector<int> R, std::vector<int> P, std::vector<int> X)
{
	/// Based on vectors Reporter, Pending, eXcluded.
	/// BronKerbosch1(R, P, X):
	///
	///   if P and X are both empty:
	///       report R as a maximal clique
	///
	///   for each vertex v in P:
	///       BronKerbosch1(R U {v}, P /\ N(v), X /\ N(v))
	///       P := P \ {v}
	///       X := X U {v}

	if(verbosity > 1) std::cout << "***************" << std::endl;

	/// if P and X are both empty
	if (P.empty() && X.empty())
	{
		/// report R as a maximal clique
		if(static_cast<int>(R.size()) >= Kmin)
			G1.maximalCliqueSet.push_back(R);
		if(verbosity > 1) std::cout << "R is maximal clique!" << std::endl;
	}

	/// Copy of P, where all operations are performed. (P is used only to get the "v"ertex)
	std::vector<int> P_copy = P;

	/// For each vertex in P
	for (std::size_t i = 0; i < P.size(); ++i)
	{
		int v = P[i];

		std::vector<int> P_intersect_Nv = {};
		std::vector<int> X_ntersect_Nv = {};

		if(verbosity > 1) 
		{
			std::cout << "(v="<<v<<") R: [";
			for(auto r:R) std::cout << r << " "; 
			std::cout << "]" << std::endl;
			std::cout << "(v="<<v<<") P: [";
			for(auto p:P) std::cout << p << " "; 
			std::cout << "]" << std::endl;
			std::cout << "(v="<<v<<") X: [";
			for(auto x:X) std::cout << x << " "; 
			std::cout << "]" << std::endl;
			std::cout << "(v="<<v<<") N: [";
			for(auto nb:neighbours[v]) std::cout << nb << " "; 
			std::cout << "]" << std::endl;
		}

		/// For each "neighbour" of v
		for (int n : neighbours[v])
		{
			/// P intersection N(v)
			for (int p : P_copy)
			{
				if (n == p)
					P_intersect_Nv.push_back(n);
			}

			/// X intersection N(v)
			for (int x : X)
			{
				if (n == x)
					X_ntersect_Nv.push_back(n);
			}
		}

		/// R union {v}
		R.push_back(v);
		/*if(verbosity > 1) {
			std::cout << "(v="<<v<<") R U {" << v << "}" << std::endl;
			std::cout << "(v="<<v<<") P_N: [";
			for(auto i1:P_intersect_Nv) std::cout << i1 << " "; 
			std::cout << "]" << std::endl;
			std::cout << "(v="<<v<<") X_N: [";
			for(auto i2:P_intersect_Nv) std::cout << i2 << " "; 
			std::cout << "]" << std::endl;
		}*/

		BronKerbosch_naive(G1, Kmin, R, P_intersect_Nv, X_ntersect_Nv);
		
		/// Set R back to previous
		R.pop_back();
		if(verbosity > 1) std::cout << "(v="<<v<<") R - {" << v << "}" << std::endl;

		/// P \ {v}
		P_copy.erase(std::remove_if(P_copy.begin(), P_copy.end(),
								[&captured_v = v](const int& x) { 
									return x == captured_v;
								}), P_copy.end());
		if(verbosity > 1) std::cout << "(v="<<v<<") P - {" << v << "}" << std::endl;

		/// X union {v}
		X.push_back(v);
		if(verbosity > 1) std::cout << "(v="<<v<<") X U {" << v << "}" << std::endl;
	}
}

void BronKerbosch_tomita(Graph &G1, int Kmin, std::vector<int> R, std::vector<int> P, std::vector<int> X)
{
	/// Based on vectors Reporter, Pending, eXcluded.
	/// BronKerbosch2(R, P, X):
	///
	///   if P and X are both empty:
	///       report R as a maximal clique
	///
	///	  choose a pivot u from P U X, to maximize |P /\ N(u)|
	///
	///   for each vertex v in P \ N(u) :
	///       BronKerbosch2(R U {v}, P /\ N(v), X /\ N(v))
	///       P := P \ {v}
	///       X := X U {v}
	
	/// if P and X are both empty
	if (P.empty() && X.empty())
	{
		/// report R as a maximal clique
		if(static_cast<int>(R.size()) >= Kmin)
			G1.maximalCliqueSet.push_back(R);
	}

	/// Choose the pivot "u" from  P U X
	int u = pivot(G1, P, X);

	/// P \ N(u) 
	std::vector<int> P_Nu;
	for (int p : P)
	{
		int check = 0;
		for (int n : neighbours[u])
			if (n == p)
				check++;				
		if(check == 0)
			P_Nu.push_back(p);
	}

	/// For each vertex in P \ N(u)
	for (std::size_t i = 0; i < P_Nu.size(); ++i)
	{
		int v = P_Nu[i];

		std::vector<int> P_intersect_Nv = {};
		std::vector<int> X_ntersect_Nv = {};

		/// For each "neighbour" of v
		for (int n : neighbours[v])
		{
			/// P intersection N(v)
			for (int p : P)
			{
				if (n == p)
					P_intersect_Nv.push_back(n);
			}

			/// X intersection N(v)
			for (int x : X)
			{
				if (n == x)
					X_ntersect_Nv.push_back(n);
			}
		}

		/// R union {v}
		R.push_back(v);

		BronKerbosch_tomita(G1, Kmin, R, P_intersect_Nv, X_ntersect_Nv);
		
		/// Set R back to previous
		R.pop_back();

		/// P \ {v}
		P.erase(std::remove_if(P.begin(), P.end(),
								[&captured_v = v](const int& x) { 
									return x == captured_v;
								}), P.end());

		/// X union {v}
		X.push_back(v);
	}
}

int pivot(Graph& G1, std::vector<int>& P, std::vector<int>& X)
{
	/// P U X
	std::vector<int> P_U_X;
	for (auto it = P.begin(); it != P.end(); ++it)
		P_U_X.push_back(*it);
	for (auto it = X.begin(); it != X.end(); ++it)
		P_U_X.push_back(*it);

	/// Pick the vertex with the lowest-degeneracy
	int picked = 0;
	for (auto it = P_U_X.begin(); it != P_U_X.end(); ++it)
		if ( picked == 0 || G1.adjList[picked].size() < G1.adjList[*it].size() )
			picked = *it;
		
	return picked;

	/// just another way
	/*int P_idx = 0;
	int P_min = G1.E;
	for (int i = P.size(); --i;)
	{
		int deg = G1.adjList[i].size();
		if ( deg < P_min )
		{
			P_min = deg;
			P_idx = i;
		}
	}

	int X_idx = 0;
	int X_min = G1.E;
	for (int i = X.size(); --i;)
	{
		int deg = G1.adjList[i].size();
		if ( deg < X_min )
		{
			X_min = deg;
			X_idx = i;
		}
	}

	return (P_min < X_min)? P_idx : X_idx;*/
}

/// O(d n 3^{d/3}) ;where [n: number of vertices of the graph] [d: degeneracy of the graph]
void BronKerbosch_eppstein(Graph &G1, int Kmin, std::vector<int> R, std::vector<int> P, std::vector<int> X)
{
	/// Degeneracy (V, E)
	/// 
	/// for each vertex vi in a degeneracy ordering v0, v1,...vn-1 of (V,E) do
	///       P := N(vi) /\ {vi+1, ...,vn-1}
	///       X := N(vi) /\ {v0,...,vi-1}
	///		  R := {vi}
	///       BronKerbosch2(R, P, X)

	/// For each vertex vi in degeneracy ordering v0,...,vn-1
	std::vector<int> degVerOrder = DegeneracyOrder(G1);
	for(std::size_t i = 0; i < degVerOrder.size(); ++i)
	{
		/// vi
		int v = degVerOrder[i];

		/// Empty for current iteration
		P.clear();
		X.clear();
		R.clear();

		/// For each "neighbour" of v
		for (std::size_t n = 0; n < neighbours[v].size(); ++n)
		{
			///  P := N(vi) /\ {vi+1, ...,vn-1}
			for (std::size_t p = i+1; p < degVerOrder.size(); ++p)
			{
				if (neighbours[v][n] == degVerOrder[p])
					P.push_back(neighbours[v][n]);
			}

			/// X := N(vi) /\ {v0,...,vi-1}
			for (std::size_t x = 0; x < i; ++x)
			{
				if (neighbours[v][n] == degVerOrder[x])
					X.push_back(neighbours[v][n]);
			}
		}

		/// R := {vi}
		R.push_back(v);

		/// BronKerbosch2(R, P, X)
		BronKerbosch_tomita(G1, Kmin, R, P, X);
	}
}

std::vector<int> DegeneracyOrder(Graph& G1)
{
	std::vector<int> deg_vertices;
	std::vector<int> deg_values;
	int j, key_values, key_vertices;

	for (int i = 0; i < G1.V; ++i)
	{
		deg_values.push_back(static_cast<int>(G1.adjList[i].size()));
		deg_vertices.push_back(i);

		key_values = deg_values[i];
		key_vertices = deg_vertices[i];
		
		j = i - 1;

		/// Shift by one position ahead the elements greater than key
		while (j >= 0 && deg_values[j] > key_values)
		{
			deg_values[j + 1] = deg_values[j];
			deg_vertices[j + 1] = deg_vertices[j];
			j = j - 1;
		}
		deg_values[j + 1] = key_values;
		deg_vertices[j + 1] = key_vertices;
	}

	/*int j, key;
	for (int i = 1; i < G1.V; ++i)
	{
		key = deg_values[i];
		j = i - 1;

		/// Move elements of arr[0..i-1],
		/// that are greater than key, to one position ahead of their current position
		while (j >= 0 && deg_values[j] > key)
		{
			deg_values[j + 1] = deg_values[j];
			j = j - 1;
		}
		deg_values[j + 1] = key;
	}*/

	return deg_vertices;
}