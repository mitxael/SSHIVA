/// ******************************************************************
///	File	: Homology.cpp
///	About	: Homology methods
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// -  https://docs.oracle.com/javase/7/docs/technotes/guides/jni/spec/functions.html
/// ******************************************************************

#include "Homology.h"


/// ************ HOMOLOGY CLASS **********************************

void Homology::PersistentHomology(Graph &G, int phEngine, int complexType, bool onlyInfinite, int dimension, bool saveCycles)
{
	switch (phEngine)
	{
		///TODO: according to complexType, calculate Filtration, then it will be required by each engine to generically-add the filtration
		///buildComplex(complexType);

		/// JAVAPLEX
		case 0:
		{
			PH_Javaplex javaplex;

			/// Load and initialize the Java Virtual Machine
			JVM JVMinstance;
			JVMinstance.StartJVM();

			/// Build Filtered-Chain-Complex from graph (0: AllCliques; 1: Graph; 2: MaxCliques, 3:MinCycles)
			//int complexType = 0;	
			jobject stream = javaplex.j_buildComplex(G, complexType, dimension);

			/// Show filtration
			if(verbosity > 0) showFiltration(G);

			/// Compute the Boundary Matrix for the complex
			/// Represented as a "Sparse Matrix" of p-chains (formal sum of p-simplices)
			/// "-1 [3] [3,5]" means that the entry between 0-simplex [3] and 1-simplex [3,5] is "-1"
			///std::map<std::vector<std::vector<int>>, double> boundaryMatrix = javaplex.createBoundaryMatrixAsDoubleSum(stream);

			/// Calculate Persistent Homology barcodes
			std::pair<jobject, jobject> PH = javaplex.computePH(G, stream);

			/// Get Betti Barcodes
			//bool onlyInfinite = true;
			bool getSequence = false;
			std::vector<int> bettiSequence = javaplex.getBettiNumbers(PH, onlyInfinite, getSequence);

			/// Get intervals (and show them)
			//int dimension = -1;		/// n: n-dimension; -1:all-dimensions
			int displayMode = 0;	/// 0: multiline 2: singleline
			//bool saveCycles = true;
			javaplex.getIntervalsAndGenerators(G, PH, dimension, displayMode, onlyInfinite, saveCycles);

			///Do not StopJVM to allow re-using
			///JVMinstance.StopJVM();
			break;	
		}

		/// Gudhi
		case 1:
		{
			///TODO: Only AllCliques works, others are not simplices!! (not complete...)
			/// Graph: expand it!	MaxCliques: Expand them!   MinCycles: also!
			buildComplex(G, complexType, dimension);

			/// Compute PH as CRWF
			PH_GUDHI::computePH(G, onlyInfinite, dimension, saveCycles);
			
			/// Set number of "connected components"
			/*int inf_at_dim0 = 0;
			for (auto& const interval : G.PH_barcodes[0])
			{
				if (interval.first.second == std::numeric_limits<double>::infinity() )
					inf_at_dim0++;
			}
			G.Components = inf_at_dim0;*/

			break;	
		}

		/// PHAT
		case 2:
		{		

			/// AllCliquesComplex(G);
			buildComplex(G, 0, dimension);

			/// Show filtration
			if(verbosity > 0) showFiltration(G);

			/// Compute PH
			PH_PHAT phat;
			phat.computePH(G, 1);
			//phat.compute_pairing(false, true, false);

			break;	
		}
		/// DiPHA
		/*case 3:
		{
			std::cout << std::endl << "********************************************" << std::endl;
			std::cout << std::endl << "**** THIS FEATURE IS UNDER CONSTRUCTION ****" << std::endl;
			std::cout << std::endl << "**** RESULTS MIGHT NOT BE FULLY CORRECT ****" << std::endl;
			std::cout << std::endl << "********************************************" << std::endl;

			/// Compute PH
			///PH_DIPHA dipha;
			///dipha.computePH(1);
			
			break;	
		}*/

		/// Perseus
		/*case 4:
		{
			std::cout << std::endl << "********************************************" << std::endl;
			std::cout << std::endl << "**** THIS FEATURE IS UNDER CONSTRUCTION ****" << std::endl;
			std::cout << std::endl << "**** RESULTS MIGHT NOT BE FULLY CORRECT ****" << std::endl;
			std::cout << std::endl << "********************************************" << std::endl;

			/// Compute PH
			///PH_PERSEUS::computePH(1);

			break;	
		}*/

		///TODO: show all PH information here, taking it from PH_barcodes and PH_filtration

	}
}

void Homology::buildComplex(Graph &G, int complexType, int dimension)
{
	/// Build EdgeReference
	ReferenceArray_ini(G, true);

	/// Build WeightRank
	bool nonDecreasing = false;
	bool zeroBased = true;
	G.buildWeightRank(nonDecreasing, zeroBased);

	/// Build simplicial complex
	switch (complexType) {
		case 0:
		{
			/// Add All Cliques to stream
			/// Persistent homology: return <MCB?>
			bool autoAlgorithm = true;
			buildCliquesComplex(G, autoAlgorithm, dimension);

			break;
		}
		case 1:					
		{
			/// Add Graph to stream
			/// Persistent homology: return <Elementary cycles> (AKA Holes AKA Hamiltonian) => minimal subgraph such that every node has degree 2
			///TODO: For G_008_complex some cycles are only expanded (e.g. 2 3 6 + 0 1; 1 5 7 + 0 2)
			buildGraphComplex(G);

			break;
		}
		case 2:	
		{
			/// Add Maximal Cliques to stream
			/// Persistent homology: return <Elementary cycles>
			///TODO: For G_008_complex & G_012_cliques results are the same (9 & 15) as buildGraphComplex()
			buildMaxCliqueComplex(G);
			
			break;
		}
		case 3:	
		{
			/// Add Minimum Cliques to stream
			/// Persistent homology: return <?>
			buildMinCyclesComplex(G);
			
			break;
		}
	}
}

void Homology::buildGraphComplex(Graph &G)
{
	/// Add vertex:		0-simplices
	for(int v = 0; v < G.V; ++v)
	{
		/// Add vertex (only javaplex)
		//...
		/// Add element
		G.PH_Filtration[0].insert(std::make_pair(std::vector<int>{v}, 0));
	}

	/// Add edges:		1-simplices
	for(int i = 0; i < G.edgesReference.size(); ++i)
	{
		int u = std::get<0>(G.edgesReference[i]);
		int v = std::get<1>(G.edgesReference[i]);
		double w = std::get<2>(G.edgesReference[i]);
		//int edgeIndex = i;	/// "i" because edges are taken in order (nondecreasing)
		int filter = G.edgeRank(u, v);
		G.PH_Filtration[filter].insert(std::make_pair(std::vector<int>{u, v}, w));
		///addElement(stream, std::vector<int>{u,v}, filter);
	}
}

void Homology::buildMaxCliqueComplex(Graph &G)
{
	/// Clique algorithm (k>=1) is chosen based on graph density
	if (G.density() > 0.4)
		BronKerbosch(G, 1, 1);	///1:Tomita (better for dense graphs)
	else
		BronKerbosch(G, 2, 1);	///2:Epstein (better for sparse graphs)

	/// Add vertices (only javaplex)
	/*for(int v = 0; v < G.V; ++v)
	{
		G.PH_Filtration[0].insert(std::make_pair(std::vector<int>{v}, 0));
	}*/

	/// Add maximal cliques
	for(auto it = G.maximalCliqueSet.begin(); it != G.maximalCliqueSet.end(); ++it)
		//for(int c = 0; c < G.cliqueS.size(); ++c)
	{
		std::vector<int> &clique = (*it);	//G.cliqueS[c];
		int k = (int)clique.size();
		if(k == 1)
		{
			/// Add 0-simplices
			int v = clique.front();
			G.PH_Filtration[0].insert(std::make_pair(std::vector<int>{v}, 0));
		}
		else if(k == 2)
		{
			/// Add 1-simplices
			int u = clique.front();
			int v = clique.back();
			G.PH_Filtration[G.edgeRank(u, v)].insert(std::make_pair(std::vector<int>{u, v}, G.edgeWeight(u,v)));
		}
		else
		{
			for(int i = 0; i < k; ++i)
			{
				int u = clique[i];
				for(int j = i + 1; j < k; ++j)
				{
					int v = clique[j];
					/// Add n-simplices
					G.PH_Filtration[G.edgeRank(u, v)].insert(std::make_pair(std::vector<int>{u, v}, G.edgeWeight(u,v)));
				}
			}
		}
	}
}

void Homology::buildMinCyclesComplex(Graph &G)
{
	/// Compute MCB
	Amaldi(G);

	/// Add vertices
	for(int v = 0; v < G.V; ++v)
	{
		G.PH_Filtration[0].insert(std::make_pair(std::vector<int>{v}, 0));
	}

	/// Add minimum cycles
	for(auto it = G.minimumCycleBasis.begin(); it != G.minimumCycleBasis.end(); ++it)
	{
		std::vector<ut> &cycle = (*it);
		int csize = (int)cycle.size();

		std::vector<int> cycleInt;
		for(auto u : cycle)
		{
			cycleInt.push_back(u);
		}
		///TODO: Add with filtration # of lowest edge?? (like WRCF)
		G.PH_Filtration[0].insert(std::make_pair(cycleInt, 0));
	}
}

void Homology::buildCliquesComplex(Graph &G, bool autoAlgorithm, int dimension)
{
	/// cliques of all dimensions
	if (dimension == -1)
		buildAllCliquesComplex(G, autoAlgorithm);
	/// only cliques up to specific dimension
	else
		allCliquesBF(G, dimension);
}

void Homology::buildAllCliquesComplex(Graph &G, bool autoAlgorithm)
{

	/// Clique algorithm (k>=1) is chosen based on graph density
	int verbold = verbosity;
	verbosity = 0;
	int kMin = 1;
	if (autoAlgorithm)
	{
		if (G.density() > 0.4)
			BronKerbosch(G, 1, kMin);	///1:Tomita (better for dense graphs)
		else
			BronKerbosch(G, 2, kMin);	///2:Eppstein (better for sparse graphs)
	}
	else
	{
		Criccaldi(G);					///3:Criccaldi
		G.cycleSpace.clear();
		G.cycleSpace_w.clear();
	}
	verbosity = verbold;

	///<K, <rank, weight, clique>[]>>
	std::map<int, std::set<std::tuple<int, double, std::vector<int>>>> discoveredCliques;

	/// Compute subcliques
	/// Each k-clique has 2^{k}-1 subcliques. e.g. k1:1  k2:3  k3:7  k4:15...
	/// For k>=1: 2^{k}-1-k subcliques. e.g. k1:0  k2:1  k3:4  k4:11...

	/// Show expected number of cliques (actual number might be lower if a clique forms part of different bigger cliques)
	/*int expectedCliques = 0;
	for ( auto it = G.cliqueS.begin(); it != G.cliqueS.end(); ++it )
	{
		int k = (int)it->size();
		expectedCliques += (int)std::pow(2,k) - 1 - k;
	}
	if(verbosity > 1) std::cout << "***** Total cliques: expected=" << expectedCliques << " & actual=" << allCliques.size() << std::endl;*/

	/// Select method to find cliques: 1=recursive; 2=iterative (faster)
	int findClique_version = 2;
	bool parallelExec = true;
	if (findClique_version == 1)
	{
		/// Method #1 (recursive)
		auto startTime = Utils::startChrono<std::chrono::high_resolution_clock::time_point>();
		for (auto it = G.maximalCliqueSet.begin(); it != G.maximalCliqueSet.end(); ++it)
		{
			int k = (int)it->size();
			std::vector<int> &clique = (*it);
			findSubCliques_v1(G, G.PH_Filtration, clique, kMin, parallelExec);
		}
		if (verbosity > 0) std::cout << "Time(FindCliques): " << Utils::stopChrono(startTime) << "s" << std::endl;

	}
	else
	{
		/// Method #2 (iterative enhanced)
		auto startTime = Utils::startChrono<std::chrono::high_resolution_clock::time_point>();
		for (auto& clique : G.maximalCliqueSet)
		{
			int k = (int)clique.size();
			int rank = (k > 1) ? G.lowestWeightIndex(clique) : 0;
			double weight = (k > 1) ? std::get<2>(G.edgesReference[rank]) : 0;
			discoveredCliques[k].insert(std::make_tuple(rank, weight, clique));
		}
		int kMax = (int)G.maximalCliqueSet.back().size();
		findSubCliques_v2(G, G.PH_Filtration, discoveredCliques, kMax, kMin, parallelExec);
		if (verbosity > 0) std::cout << "Time(FindCliques): " << Utils::stopChrono(startTime) << "s" << std::endl;
	}

	///TODO: Filtration is added to graph by each PhEngine

}

int Homology::findSubCliques_v1(Graph &G, std::map<int, std::set<std::pair<std::vector<int>, double>>> &allCliques, 
				std::vector<int> &clique, int kMin,  bool parallelExec)
 {
	/// Get information on current clique
	int k = (int)clique.size();

	if(k < 2)
	{
		/// Store the current clique "k=1"
		int rank = 0;
		double w = 0;
		allCliques[rank].insert(std::make_pair(clique, w));
	}
	else /// (k>=2)
	{
		/// Store the current clique "k"
		int rank = G.lowestWeightIndex(clique);
		double w = std::get<2>(G.edgesReference[rank]);	
		allCliques[rank].insert(std::make_pair(clique, w));
		
		/// Recursive finding of all cliques "k-1"
		int i, sub;		

		/// OMP: make parallel
		int max_threads = omp_get_max_threads();
		#pragma omp parallel num_threads(max_threads) if (parallelExec)
		for(i = 0; i < clique.size(); ++i)
		{
			//std::cout << "Thread number: " << omp_get_thread_num() << std::endl;
			std::vector<int> subclique = clique;
			subclique.erase(subclique.begin()+i);

			/// OMP: region to be executed by a single thread
			#pragma omp single
			{
				///Regular execution
				if(subclique.size() >= kMin)
				{
					sub = findSubCliques_v1(G, allCliques, subclique, kMin, parallelExec);
				}
			}
		}
	}

	return k;
}

void Homology::findSubCliques_v2(Graph &G, std::map<int, std::set<std::pair<std::vector<int>, double>>> &allCliques, 
				std::map<int, std::set<std::tuple<int, double, std::vector<int>>>> &discoveringCliques, int kMax, int kMin, bool parallelExec)
 {
	/// TODO: Try to use only one container rather than two (discoveringCliques and allCliques)
	/// TODO: [key]: rank vs k-size

	/// This method loops each clique "k" to get all "k-1" cliques (by removing one different vertex at each iteration)
	/// Subcliques are stored in a map, where the key is the "k" number of cliques

	/// get # of available threads
	int max_threads = omp_get_max_threads();

	/// set all available threads for being used for parallelism
	//omp_set_num_threads(omax_threads);

	/// 0: threads are assigned only to the top-level parallel region (n) 1: threads are assigned at each level (n-1*i)
	//omp_set_nested(0);

	/// For each "k-set" of cliques (starting from the biggest one)
	for (int k = kMax; k > kMin; k--)
	{
		/// Make available an iterator to the first "clique_tuple" in the "k" set
		std::set<std::tuple<int, double, std::vector<int>>>::iterator it = discoveringCliques[k].begin();

			/// For each "clique_tuple" in the "k-set" (starting from the first one)
			int n_kcliques = discoveringCliques[k].size();
			for(int s = 0; s < n_kcliques; ++s)
			{	
				/// Dereference the "clique_tuples" and Extract individual data from it
				std::tuple<int, double, std::vector<int>> clique_tuple = *it;
				std::vector<int> clique = std::get<2>(clique_tuple);
				double weight = std::get<1>(clique_tuple);
				int rank = std::get<0>(clique_tuple);
				
				///Add non-vertices to ALLCLIQUES
				allCliques[rank].insert(std::make_pair(clique, weight));

				/// Get lowestEdge vertices
				/*int src = G.WeightRank[rank].first;
				int dst = G.WeightRank[rank].second;
				if(vertex!=src || vertex!=dst) new = current*/
			
				/// OMP: Distribute the processing of clique's vertices" among a team of threads
				int n_clique = clique.size();
				#pragma omp parallel for num_threads(max_threads) if(parallelExec)
				for (int v = 0; v < n_clique; ++v)
				{
					/// Create new clique as the copy of source clique BUT excluding the current "vertex"
					int& vertex = clique[v];
					std::vector<int> new_clique;
					std::copy_if(clique.begin(), clique.end(), std::back_inserter(new_clique), [vertex](const int& elem) { return elem != vertex; });				

					/// Check if removed vertex is related to lower edge
					int new_rank = (k > 1) ? G.lowestWeightIndex(new_clique) : 0;
					double new_weight = (k > 1) ? std::get<2>(G.edgesReference[new_rank]) : 0;
					int kNew = k - 1;
					
					/// OMP: Access allowed for one thread at time
					#pragma omp critical
					{
						/// Add new clique
						discoveringCliques[kNew].insert(std::make_tuple(new_rank, new_weight, new_clique));

						///Add vertices to ALLCLIQUES
						if (kNew == 1) allCliques[new_rank].insert(std::make_pair(new_clique, new_weight));
					}
					
				}

				/// Once the "clique_tuple" has been processed, move to the iterator to the next one
				it++;
			}
	}
	///bool end = true;
}

void Homology::showFiltration(Graph &G)
{
	/// Show filtration simplices
	std::cout << "Filtrations:" << std::endl;
	int counter = 0;
	int max = 0;
	for (auto it = G.PH_Filtration.begin(); it != G.PH_Filtration.end(); ++it)
	{
		const int &rank = it->first;
		const auto &rankCliques = (it->second);
		std::cout << "\tFilter<" << rank << ">: ";
		bool firstPair = true;
		for (auto cliqueP : rankCliques)
		{
			std::vector<int> &clique = cliqueP.first;

			int k = (int)clique.size();
			if(k > max) max = k;
				if (k == 1)
				{
					if (!firstPair) std::cout << " ; ";
						std::cout << "[" << clique[0] << "]";
				}
				else
				{
						if (!firstPair) std::cout << " ; ";
						std::cout << "[";
						for(auto &vertex : clique)
						{
							std::cout << vertex;
							if(vertex != clique.back())
								std::cout << ",";
						}
						std::cout << "]";
				}
			firstPair = false;
			counter++;
		}
		std::cout << std::endl;
	}

	/// Show filtration info
	std::cout << "\t*N. of vertices = " << G.V << std::endl;
	std::cout << "\t*N. of simplices = " << counter << std::endl;
	std::cout << "\t*Dimension = " << max << std::endl;
	///std::cout << "\t*Dimension = " << G.maximalCliqueSet[G.maximalCliqueSet.size()-1].size() << std::endl;
}

double Homology::bottleneckDistance(std::string file1, std::string file2, double tolerance)
{
	/// CHECK FILE#1
	const char *file1_ = file1.c_str();
	std::ifstream fin1(file1_);
	if (fin1.is_open())
	{
		/// CHECK FORMAT
		std::string values;
		while ( values.empty() && !fin1.eof() )
			std::getline(fin1, values);
		/*std::stringstream sstream(values);
		std::vector<double> line;
		while (1)
		{ 
			double n; sstream >> n; if (!sstream) break;
			line.push_back(n);
		}
		if (line.size() < 2) return -1;*/
		double line[4];
		int n = sscanf_s(values.c_str(), "%lf %lf %lf %lf", &line[0], &line[1], &line[2], &line[3]);
		if (n < 2 || n > 3) return -1;
	}
	else return -1;

	/// CHECK FILE#2
	const char *file2_ = file2.c_str();
	std::ifstream fin2(file2_);
	if (fin2.is_open())
	{
		/// CHECK FORMAT
		std::string values;
		while ( values.empty() && !fin2.eof() )
			std::getline(fin2, values);
		/*std::stringstream sstream(values);
		std::vector<double> line;
		while (1)
		{ 
			double n; sstream >> n; if (!sstream) break;
			line.push_back(n);
		}
		if (line.size() < 2) return -2;*/
		double line[4];
		int n = sscanf_s(values.c_str(), "%lf %lf %lf %lf", &line[0], &line[1], &line[2], &line[3]);
		if (n < 2 || n > 3) return -1;
	}
	else return -2;
	
	double b = PH_GUDHI::bottleneckDistance(file1, file2, tolerance);

	return b;
}