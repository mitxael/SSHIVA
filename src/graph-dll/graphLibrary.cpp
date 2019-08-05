/// ***********************************************************
///	File	: graphLibrary.cpp
///	About	: C++/CLI interface between C++ library and WPF C# GUI
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
/// ******************************************************************

#include "graphLibrary.h"
#include <msclr\marshal_cppstd.h>
//#include <cliext/list>			/// Many extended classes: lists, etc.

namespace GRAPH
{
	/// Namespaces
	using namespace System::Collections::Generic;
	//using namespace System;

	/// ************* CORE METHODS ************************************************************
	GraphManager::GraphManager(int V)
	{
		m_Graph = new Graph(V);
		getCore();
	}

	void GraphManager::buildGraph(bool save2file, int n)
	{
		/// Create unmanaged graph
		m_Graph->buildGraph(save2file, n);
		
		/// Update managed graph
		set_mGraph(*this->m_Graph);
	}

	void GraphManager::buildGraph(bool save2file,int n, int m)
	{
		/// Create unmanaged graph
		m_Graph->buildGraph(save2file, n, m);
		
		/// Update managed graph
		set_mGraph(*this->m_Graph);
	}

	void GraphManager::buildGraph(bool save2file, int n, double density)
	{
		/// Create unmanaged graph
		m_Graph->buildGraph(save2file, n, density);
		
		/// Update managed graph
		set_mGraph(*this->m_Graph);
	}

	void GraphManager::buildGraph(System::String ^filename)
	{	
		/// Parse "file" to "u_file"
		msclr::interop::marshal_context ctx;
		std::string u_file = ctx.marshal_as<std::string>(filename);
		
		/// Create unmanaged graph
		m_Graph->buildGraph(u_file);

		/// Update managed graph
		set_mGraph(*this->m_Graph);
	}

	void GraphManager::saveGraph(System::String ^filename, char format)
	{
		msclr::interop::marshal_context ctx;
		
		/// Parse "file" to "u_file"
		std::string u_file = ctx.marshal_as<std::string>(filename);

		/// Return function
		m_Graph->saveGraph(u_file, format);
	}

	void GraphManager::getCore()
	{
		this->V = m_Graph->V;
		this->E = m_Graph->E;
		this->K = m_Graph->K;
		this->Density = m_Graph->Density;
		this->Degeneracy = m_Graph->Degeneracy;
		this->Cyclomatic = m_Graph->Cyclomatic;
		this->Components = m_Graph->Components;
		this->source = gcnew System::String(m_Graph->source.c_str());
	}

	void GraphManager::getAdjList()
	{
		this->adjList = gcnew cli::array<cli::array<System::Tuple<int,double>^>^> (m_Graph->V);
		for (int i = 0; i < m_Graph->V; ++i)
		{
			this->adjList[i] = gcnew cli::array<System::Tuple<int,double>^>((int)m_Graph->adjList[i].size());
			int j = 0;
			for (auto it = m_Graph->adjList[i].begin(); it != m_Graph->adjList[i].end(); ++it)	
			{
				this->adjList[i][j] = gcnew System::Tuple<int, double>((*it).first, (*it).second);
				j++;
			}
		}
	}

	void GraphManager::getAlgorithmS()
	{
		int nAlgoS = m_Graph->algorithmS.size();
		this->algorithmS = gcnew cli::array<System::String^>(nAlgoS);
		
		msclr::interop::marshal_context ctx;

		for (int i = 0; i < nAlgoS; ++i)
		{
			/// Parse "u_string" to "m_string"
			this->algorithmS[i] =  gcnew System::String(m_Graph->algorithmS[i].c_str());
		}
	}

	void GraphManager::getEdgesReference(Graph &u_G1)
	{
		this->EdgesReference = gcnew cli::array<System::Tuple<int,int, double>^>((int)u_G1.edgesReference.size());
		for (int i = 0; i < u_G1.edgesReference.size(); ++i)
		{
			int src = std::get<0>(u_G1.edgesReference[i]);
			int dst = std::get<1>(u_G1.edgesReference[i]);
			double wgt = std::get<2>(u_G1.edgesReference[i]);
			this->EdgesReference[i] = gcnew System::Tuple<int, int, double>(src, dst, wgt);
		}
	}

	void GraphManager::getCycleSpace(Graph &u_G1)
	{
		this->cycleSpace = gcnew cli::array<cli::array<int>^>((int)u_G1.cycleSpace.size());
		for (int i = 0; i < u_G1.cycleSpace.size(); ++i)
		{
			this->cycleSpace[i] = gcnew cli::array<int>((int)u_G1.cycleSpace[i].size());
			for (int j = 0; j < u_G1.cycleSpace[i].size(); ++j)
			{
				this->cycleSpace[i][j] = u_G1.cycleSpace[i][j];
			}
		}

		this->minimumCycleBasis_w = gcnew cli::array<double>((int)u_G1.minimumCycleBasis_w.size());
		for (int k = 0; k < u_G1.minimumCycleBasis_w.size(); ++k)
		{
				this->minimumCycleBasis_w[k] = u_G1.minimumCycleBasis_w[k];
		}
	}

	void GraphManager::getMinimumCycles(Graph &u_G1)
	{
		this->minimumCycleBasis = gcnew cli::array<cli::array<int>^>((int)u_G1.minimumCycleBasis.size());
		for (int i = 0; i < u_G1.minimumCycleBasis.size(); ++i)
		{
			this->minimumCycleBasis[i] = gcnew cli::array<int>((int)u_G1.minimumCycleBasis[i].size());
			for (int j = 0; j < u_G1.minimumCycleBasis[i].size(); ++j)
			{
				this->minimumCycleBasis[i][j] = u_G1.minimumCycleBasis[i][j];
			}
		}

		this->minimumCycleBasis_w = gcnew cli::array<double>((int)u_G1.minimumCycleBasis_w.size());
		for (int k = 0; k < u_G1.minimumCycleBasis_w.size(); ++k)
		{
				this->minimumCycleBasis_w[k] = u_G1.minimumCycleBasis_w[k];
		}
	}

	void GraphManager::getMaximalCliques(Graph &u_G1)
	{
		this->maximalCliqueSet = gcnew cli::array<cli::array<int>^>((int)u_G1.maximalCliqueSet.size());
		for (int i = 0; i < u_G1.maximalCliqueSet.size(); ++i)
		{
			this->maximalCliqueSet[i] = gcnew cli::array<int>((int)u_G1.maximalCliqueSet[i].size());
			for (int j = 0; j < u_G1.maximalCliqueSet[i].size(); ++j)
			{
				this->maximalCliqueSet[i][j] = u_G1.maximalCliqueSet[i][j];
			}
		}
	}

	void GraphManager::getCliqueSpace(Graph &u_G1)
	{
		/// Define cliqueSpace
		int spaceSize = (int)u_G1.cliqueSpace.size();
		this->cliqueSpace = gcnew Dictionary<int, List<cli::array<int>^>^>();

		/// For each rank in the Map
		for (auto const& mapEntry : u_G1.cliqueSpace)
		{
			/// Get rank
			int m_rank = mapEntry.first;

			/// Get cliques
			std::set<std::vector<int>> u_cliqueS = mapEntry.second;
			int cliquesSize = (int)u_cliqueS.size();
			
			System::Collections::Generic::List<cli::array<int>^>^ m_cliqueS = 
										gcnew System::Collections::Generic::List<cli::array<int>^>();
			
			/// For each Clique in the Set
			for (auto const& u_clique : u_cliqueS)
			{
				/// Clique
				int cliqueSize = (int)u_clique.size();
				cli::array<int>^ m_clique = gcnew cli::array<int>(cliqueSize);
				for (int i = 0; i < cliqueSize; ++i)
				{
					m_clique[i] = u_clique[i];
				}

				/// Add clique
				m_cliqueS->Add(m_clique);

			}
			/// Add current dimension
			this->cliqueSpace->Add(m_rank, m_cliqueS);
		}
	}

	void GraphManager::getPhFiltration(Graph &u_G1)
	{
		/// Define PH_Filtration
		int spaceSize = (int)u_G1.PH_Filtration.size();
		this->PH_Filtration = gcnew Dictionary<int, List<System::Tuple<cli::array<int>^, double>^>^>();

		/// For each rank in the Map
		for (auto const& mapEntry : u_G1.PH_Filtration)
		{
			/// Get rank
			int m_rank = mapEntry.first;

			/// Get filtrations
			std::set<std::pair<std::vector<int>, double>> u_filtrationS = mapEntry.second;
			int filtrationsSize = (int)u_filtrationS.size();
			
			System::Collections::Generic::List<System::Tuple<cli::array<int>^, double>^>^ m_filtrationS = 
										gcnew System::Collections::Generic::List<System::Tuple<cli::array<int>^, double>^>();
			
			/// For each Clique in the Set
			for (auto const& u_filtrationTuple : u_filtrationS)
			{
				/// Clique
				std::vector<int> u_filtration = u_filtrationTuple.first;
				int filtrationSize = (int)u_filtration.size();
				cli::array<int>^ m_filtration = gcnew cli::array<int>(filtrationSize);
				for (int i = 0; i < filtrationSize; ++i)
				{
					m_filtration[i] = u_filtration[i];
				}

				/// Rank
				double m_weight = u_filtrationTuple.second;	

				/// Add filtration tuple
				System::Tuple<cli::array<int>^, double>^ m_filtrationTuple = gcnew System::Tuple<cli::array<int>^, double>(m_filtration, m_weight);
				m_filtrationS->Add(m_filtrationTuple);

			}
			/// Add current dimension
			this->PH_Filtration->Add(m_rank, m_filtrationS);
		}
	}

	void GraphManager::getPhBarcodes(Graph &u_G1)
	{
		/// Define m_PH_barcodes
		int phSize = (int)u_G1.PH_barcodes.size();
		///this->PH_barcodes = gcnew Dictionary<int,Dictionary<System::Tuple<double,double>^, List<cli::array<cli::array<int>^>^>^>^>();
		this->PH_barcodes = gcnew Dictionary< int, System::Collections::Generic::List<System::Collections::Generic::KeyValuePair<
						System::Tuple<double, double>^,
						System::Collections::Generic::List<cli::array<cli::array<int>^>^>^
						>>^>();

		/// For each dimension in u_PH_barcodes
		for (auto const& barcodesAtDimensionN : u_G1.PH_barcodes)
		{
			/// Get dimension
			int dim = barcodesAtDimensionN.first;

			/// Get barcodesCollection
			std::multimap<std::pair<double,double>, std::vector<std::vector<int>>> collection = barcodesAtDimensionN.second;
			int collSize = (int)collection.size();
			///Dictionary<System::Tuple<double, double>^, List<cli::array<cli::array<int>^>^>^>^ m_collection = 
			///							gcnew Dictionary<System::Tuple<double, double>^, List<cli::array<cli::array<int>^>^>^>();
			List<KeyValuePair<
						System::Tuple<double, double>^,
						List<cli::array<cli::array<int>^>^>^
						>>^ m_collection = gcnew List<KeyValuePair<
						System::Tuple<double, double>^,
						List<cli::array<cli::array<int>^>^>^
						>>();

			System::Tuple<double, double>^ last_interval = gcnew System::Tuple<double, double>(-1, -1);
			for (auto const& barcode : collection)
			{
				/// Interval
				std::pair<double, double> u_interval = barcode.first;
				System::Tuple<double, double>^ m_interval = gcnew System::Tuple<double, double>(u_interval.first, u_interval.second);
				
				/// Generators
				std::vector<std::vector<int>> u_generators = barcode.second;
				int genSize = (int)u_generators.size();
				cli::array<cli::array<int>^>^ m_generators = gcnew cli::array<cli::array<int>^>(genSize);
				for (int i = 0; i < genSize; ++i)
				{
					std::vector<int> &u_element = u_generators[i];
					int elemSize = (int)u_element.size();
					cli::array<int>^ m_element = gcnew cli::array<int>(elemSize);
					for (int j = 0; j < elemSize; ++j)
					{
						m_element[j] = u_element[j];
					}
					m_generators[i] = m_element;
				}

				/// Check if repeated interval
				List<cli::array<cli::array<int>^>^>^ list_generators = gcnew List<cli::array<cli::array<int>^>^>;

				/// Commented because it's not more required to prevent duplicates: List<KeyValuePair>
				/*if(m_interval->Item1 != last_interval->Item1 || m_interval->Item2 != last_interval->Item2)
				{*/
					/// Add generators to the List of current interval
					list_generators->Add(m_generators);

					/// Add to collection of current dimension
					///m_collection->Add(m_interval, list_generators);
					KeyValuePair<System::Tuple<double, double>^, List<cli::array<cli::array<int>^>^>^> kp(m_interval, list_generators);
					m_collection->Add(kp);
				/*}
				else
				{
					/// Add generators to existing List of last_interval
					List<cli::array<cli::array<int>^>^>^ last_list = gcnew List<cli::array<cli::array<int>^>^>;
					///m_collection->TryGetValue(m_interval, last_list);	/// get element with "key==m_interval" into last_list
					for each (KeyValuePair<System::Tuple<double, double>^,List<cli::array<cli::array<int>^>^>^> pair in m_collection)
					{
						if(pair.Key->Item1 == m_interval->Item1 && pair.Key->Item2 == m_interval->Item2)
							last_list = pair.Value;
					}
					last_list->Add(m_generators);
				}

				///Store last interval
				last_interval = m_interval;*/
			}
			/// Add current dimension
			this->PH_barcodes->Add(dim, m_collection);
		}
	}

	void GraphManager::reset()
	{
		return m_Graph->reset();
	}

	void GraphManager::set_mGraph(Graph& u_Graph)
	{	
		/// Get this->*
		this->getCore();
		this->getAdjList();
		this->getAlgorithmS();
		this->getCycleSpace(u_Graph);
		this->getCliqueSpace(u_Graph);
		this->getEdgesReference(u_Graph);
		this->getMaximalCliques(u_Graph);
		this->getMinimumCycles(u_Graph);
		this->getPhBarcodes(u_Graph);
		this->getPhFiltration(u_Graph);

		/// Get this->m_graph
		*(this->m_Graph) = u_Graph;
	}

	Graph GraphManager::set_uGraph()
	{
		Graph u_Graph = Graph(this->m_Graph->V);
		u_Graph = (*this->m_Graph);
		/*u_Graph.adjList = this->m_Graph->adjList;
		u_Graph.algorithmS = this->m_Graph->algorithmS;
		u_Graph.BFS_spanningTree = this->m_Graph->BFS_spanningTree;
		u_Graph.cliqueSpace = this->m_Graph->cliqueSpace;
		u_Graph.Components = this->m_Graph->Components;
		u_Graph.cycleSpace = this->m_Graph->cycleSpace;
		u_Graph.cycleSpace_w = this->m_Graph->cycleSpace_w;
		u_Graph.Cyclomatic = this->m_Graph->Cyclomatic;
		u_Graph.Density = this->m_Graph->Density;
		u_Graph.Degeneracy = this->m_Graph->Degeneracy;
		u_Graph.DFS_spanningTree = this->m_Graph->DFS_spanningTree;
		u_Graph.E = this->m_Graph->E;
		u_Graph.edgesReference = this->m_Graph->edgesReference;
		u_Graph.K = (u_Graph.K != 3)? this->K : this->m_Graph->K;	/// Get K as: <parameter> OR <default=3>
		u_Graph.maximalCliqueSet = this->m_Graph->maximalCliqueSet;
		u_Graph.minimumCycleBasis = this->m_Graph->minimumCycleBasis;
		u_Graph.minimumCycleBasis_w = this->m_Graph->minimumCycleBasis_w;
		u_Graph.PH_barcodes = this->m_Graph->PH_barcodes;
		u_Graph.PH_Filtration = this->m_Graph->PH_Filtration;
		u_Graph.source = this->m_Graph->source;
		u_Graph.V = this->m_Graph->V;
		u_Graph.vTags = this->m_Graph->vTags;
		u_Graph.weightRank = this->m_Graph->weightRank;*/
		
		return u_Graph;
	}

	/// ************** UTILITY METHODS **********************************************************

	double UtilsManager::measurePerformance(GraphManager ^%G1, int algorithm, int verbosity, bool hiresClock, bool saveToFile, listOfParams algoParams)
	{
		/// Prepare for unmanaged processing
		Graph u_G1 = G1->set_uGraph();

		/// Execute algorithm
		double timing = m_Utils->measurePerformance(u_G1, algorithm, verbosity, hiresClock, saveToFile, convert_mParams2uParams(algoParams));

		/// Prepare for managed handling
		G1->set_mGraph(u_G1);
		G1->verbosity = verbosity;

		return timing;
	}

	std::vector<std::pair<std::string, std::string>> UtilsManager::convert_mParams2uParams(listOfParams algoParams)
	{
		std::vector<std::pair<std::string, std::string>> u_algoParams;
		msclr::interop::marshal_context ctx;
		for each (System::Tuple<System::String^, System::String^>^ param in algoParams)
		{
			/// Parse "file" to "u_file"
			std::string u_param_key = ctx.marshal_as<std::string>(param->Item1);
			std::string u_param_value = ctx.marshal_as<std::string>(param->Item2);
			std::pair<std::string, std::string> pair = std::make_pair(u_param_key, u_param_value);
			u_algoParams.push_back(pair);
		}	

		return u_algoParams;
	}

	int UtilsManager::outputToFile_Take()
	{
		return m_Utils->outputToFile_Take();
	}

	int UtilsManager::outputToFile_Release(int old_stdout)
	{
		return m_Utils->outputToFile_Release(old_stdout);
	}

	std::string UtilsManager::setLogFilename()
	{
		return m_Utils->setLogFilename();
	}

	void UtilsManager::terminateExecution()
	{
		return m_Utils->terminateExecution();
	}

	std::string UtilsManager::convert_mString2uString(System::String^ m_string)
	{
		msclr::interop::marshal_context ctx;
		std::string u_string = ctx.marshal_as<std::string>(m_string);

		return u_string;
	}

	char UtilsManager::convert_mChar2uChar(System::Char^ m_char)
	{
		System::String^ m_string = "";
		m_string += m_char;

		msclr::interop::marshal_context ctx;
		wchar_t * u_string = (wchar_t *)ctx.marshal_as<wchar_t const *>(m_string);
		char u_char = (char)u_string[0];

		return u_char;
	}

	/// ************** HOMOLOGY METHODS **********************************************************
	
	///void Homology::bottleneckDistance(std::string file1, std::string file2, double tolerance)
	double HomologyManager::bottleneckDistance(System::String^ file1, System::String^ file2, double tolerance)
	{	
		/// Parse "files" to "u_files"
		msclr::interop::marshal_context ctx;
		std::string u_file1 = ctx.marshal_as<std::string>(file1);
		std::string u_file2 = ctx.marshal_as<std::string>(file2);

		if(tolerance <= 0)
			tolerance = std::numeric_limits<double>::min();

		/// Execute algorithm
		double b = m_Homology->bottleneckDistance(u_file1, u_file2, tolerance);

		return b;
	}
}