/// ***********************************************************
///	File	: graphLibrary.h
///	About	: Header of graphLibrary class
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
/// ******************************************************************

#pragma once
#include "../graph-lib/base/Graph.h"
#include "../graph-lib/base/Utils.h"
//#include "../graph-lib/base/LinearIndependence.h"
//#include "../graph-lib/cycles/Amaldi.h"
//#include "../graph-lib/cycles/DePina.h"
//#include "../graph-lib/cycles/Horton.h"
//#include "../graph-lib/cycles/Hybrid.h"
//#include "../graph-lib/cliques/BronKerbosch.h"
//#include "../graph-lib/cliques/Criccaldi.h"
//#include "../graph-lib/cliques/Dharwadker.h"
#include "../graph-lib/homology/Homology.h"
//#include <cliext/map>

namespace GRAPH
{
	public ref class GraphManager
   {
	/// Internal assembly scope (i.e. only accessible from code in the same .exe or .dll)
	internal:
		Graph * m_Graph;				/// managed Graph

	/// Private class scope (i.e. accessible only from code in the same class)
	private:
		/// Custom types
		/*typedef System::Collections::Generic::Dictionary<
						int,
						System::Collections::Generic::Dictionary<
							System::Tuple<double, double>^,System::Collections::Generic::List<
								cli::array<cli::array<int>^>^
							>^
						>^
					>^ m_multimap;*/
		typedef System::Collections::Generic::Dictionary< 
							int, 
							System::Collections::Generic::List<System::Collections::Generic::KeyValuePair<
								System::Tuple<double, double>^,
								System::Collections::Generic::List<cli::array<cli::array<int>^>^>^
							>>^
						>^ m_multimap;
		typedef System::Collections::Generic::Dictionary<
							int,
							System::Collections::Generic::List<
								System::Tuple<cli::array<int>^, double>^
							>^
						>^ m_setmap;

	/// Public class scope (i.e. accessible from code outside the class)
	public:	
		/// members
		int V;    								
		int E;
		int K;
		double Density;
		int Degeneracy;
		int Cyclomatic;
		int Components;
		int verbosity;
		
		cli::array<cli::array<System::Tuple<int, double>^>^>^ adjList;
		cli::array<System::Tuple<int, int, double>^>^ EdgesReference;
		
		System::String^ source;
		cli::array<System::String^>^ algorithmS;

		cli::array<cli::array<int>^>^ cycleSpace;
		cli::array<cli::array<int>^>^ minimumCycleBasis;
		cli::array<double>^ minimumCycleBasis_w;

		System::Collections::Generic::Dictionary<int,System::Collections::Generic::List<cli::array<int>^>^>^ cliqueSpace;
		cli::array<cli::array<int>^>^ maximalCliqueSet;
		
		m_setmap PH_Filtration;
		m_multimap PH_barcodes;

		/// methods
		GraphManager(int V);						
		void buildGraph(bool save2file,int n);
		void buildGraph(bool save2file,int n, int m);
		void buildGraph(bool save2file,int n, double density);
		void buildGraph(System::String ^filename);
		void saveGraph(System::String ^filename, char format);
		void getCore();
		void getAdjList();
		void getAlgorithmS();
		void getEdgesReference(Graph &u_G1);
		void getCycleSpace(Graph &u_G1);
		void getMinimumCycles(Graph &u_G1);
		void getCliqueSpace(Graph &u_G1);
		void getMaximalCliques(Graph &u_G1);
		void getPhFiltration(Graph &u_G1);
		void getPhBarcodes(Graph &u_G1);
		void reset();
		void set_mGraph(Graph& u_Graph);
		Graph set_uGraph();
   };

   public ref class UtilsManager
   {
	internal:
	   Utils* m_Utils;					/// managed Utils

	private:
		/// Custom types
	   typedef System::Collections::Generic::List<System::Tuple<System::String^, System::String^>^>^ listOfParams;

   public:
		double measurePerformance(GRAPH::GraphManager ^%G1, int algorithm, int verbosity, bool hiresClock, bool saveToFile, listOfParams algoParams);
		int outputToFile_Take();
		int outputToFile_Release(int old_stdout);
		std::string setLogFilename();
		void terminateExecution();

	   	/// Convertors
		std::vector<std::pair<std::string, std::string>> convert_mParams2uParams(listOfParams algoParams);
		std::string convert_mString2uString(System::String^ m_string);
		char convert_mChar2uChar(System::Char^ m_char);

   };

	public ref class HomologyManager
   {
	internal:
	   Homology* m_Homology;					/// managed homology

	private:
		/// Custom types

   public:
		double bottleneckDistance(System::String^ file1, System::String^ file2, double tolerance);

   };
}