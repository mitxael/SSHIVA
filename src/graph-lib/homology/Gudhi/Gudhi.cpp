
#include "Gudhi.h"

#include "include/gudhi/read_persistence_from_file.h"
#include "include/gudhi/reader_utils.h"
#include "include/gudhi/graph_simplicial_complex.h"


typedef int Vertex_handle;
typedef double Filtration_value;
typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS,
		boost::property <Gudhi::vertex_filtration_t, Filtration_value >,
		boost::property <Gudhi::edge_filtration_t, Filtration_value > > Graph_t;

void PH_GUDHI::computePH(Graph& G, bool onlyInfinite, int dimension, bool saveCycles)
{
	/// START CLOCK
	auto startTime = Utils::startChrono<std::chrono::high_resolution_clock::time_point>();

	Gudhi::Simplex_tree<> complex;

	/// ALLCLIQUES by mitXael
	complexFromFiltration(G, complex);
	/// ALLCLIQUES by Gudhi
	//buildCliqueComplex(G, complex);

	/// Show filtration
	if(verbosity > 0) Homology::showFiltration(G);
	

	/// SHOW VERTICES IN COMPLEX
	/*std::cout << "Iterator on vertices: ";
	for (auto vertex : st.complex_vertex_range()) {
		std::cout << vertex << " ";
	}
	std::cout << std::endl;*/

	/// SHOW SIMPLICES IN COMPLEX
	/*std::cout << "Iterator on simplices: " << std::endl;
	for (auto simplex : st.complex_simplex_range()) {
		std::cout << "   ";
		for (auto vertex : st.simplex_vertex_range(simplex)) {
			std::cout << vertex << " ";
		}
		std::cout << std::endl;
	}*/

	/// SHOW SIMPLICES PLUS FILTRATION
	/*if(verbosity > 0)
	{
		std::cout << "Filtration:" << std::endl;
		for (auto f_simplex : complex.filtration_simplex_range()) {
			std::cout << "\t[" << complex.filtration(f_simplex) << "] ";
			for (auto vertex : complex.simplex_vertex_range(f_simplex)) {
				std::cout << vertex << " ";
			}
			std::cout << std::endl;
		}
		std::cout << "\t*N. of vertices = " << complex.num_vertices() << std::endl;
		std::cout << "\t*N. of simplices = " << complex.num_simplices() << std::endl;
		std::cout << "\t*Dimension = " << complex.dimension() << std::endl;
	}*/

	/// SHOW SIMPLICES PLUS BOUNDARIES
	/*std::cout << "Iterator on Simplices in the filtration, and their boundary simplices:" << std::endl;
	for (auto f_simplex : st.filtration_simplex_range()) {
		/// FILTRATION
		std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
		/// VERTICES
		for (auto vertex : st.simplex_vertex_range(f_simplex)) {
			std::cout << vertex << " ";
		}
		std::cout << std::endl;
		/// BOUNDARIES
		for (auto b_simplex : st.boundary_simplex_range(f_simplex)) {
			std::cout << "      " << "[" << st.filtration(b_simplex) << "] ";
			for (auto vertex : st.simplex_vertex_range(b_simplex)) {
				std::cout << vertex << " ";
			}
			std::cout << std::endl;
		}
	}*/

	/// COMPUTE PERSISTENT HOMOLOGY
	int coeff_field_characteristic = 2;		///Characteristic p of the coefficient field Z/pZ for computing homology.
	Filtration_value min_persistence = 0.0;	///Minimal lifetime of homology feature to be recorded. 
	Persistent_cohomology pcoh(complex);
	pcoh.init_coefficients(coeff_field_characteristic);
	pcoh.compute_persistent_cohomology(min_persistence);
	
	/// BETTI NUMBERS (ALL)
	std::vector<int> betti = pcoh.betti_numbers();
	std::cout << "Betti numbers: {";
	for (int i = 0; i < betti.size(); ++i)
	{
		std::cout << i << ": " << betti[i];
		if(i < betti.size()-1)
			std::cout << ", ";
	}
	std::cout << "}" << std::endl;

	/// BETTI NUMBERS (PERSISTENT)
	/*auto pbetti = pcoh.persistent_betti_numbers(0, 99);
	std::cout << "(Persistent) Betti numbers:" << std::endl;
	for(int i = 0; i < pbetti.size(); ++i)
		std::cout << "dim-" << i << ":" << pbetti[i] << std::endl;*/

	/// INTERVALS
	///Persistence_intervals phinter(vector_of_pairs);
	Barcodes intervals = pcoh.get_persistent_pairs();
	std::cout << "Intervals:" << std::endl;
	for (int i = 0; i < intervals.size(); ++i)
	{
		if( (dimension == -1) || (dimension == i) )
		{
			auto intervals_at = pcoh.intervals_in_dimension(i);
			std::sort(intervals_at.begin(), intervals_at.end(), [](auto &left, auto &right) {
				return left.first < right.first;
			});
			if (intervals_at.size() > 0)
			{
				std::cout << "\tdim-" << i << ":" << std::endl;
				for (auto&& pair : intervals_at)
				{
					if( (!onlyInfinite) || (onlyInfinite && (pair.second == std::numeric_limits<double>::infinity())) )
						std::cout << "\t[" << pair.first << " , " << pair.second << ") " << std::endl;
				}
			}
		}
	}

	/// DIAGRAM
	/*std::cout << "Z : dim : birth : death" << std::endl;
	pcoh.output_diagram();*/

	/// Add barcodes to graph 
	barcodesIntoGraph(G, pcoh);

	/// Add cycles to graph (NOT SUPPORTED BECAUSE THERE ARE NO GENERATORS)
	///if(saveCycles)
	///...

	/// END CLOCK
	if(verbosity>0) std::cout << "Time(computePH): " << Utils::stopChrono(startTime) << "s" << std::endl;
}

void PH_GUDHI::complexFromFiltration(Graph& G,  Gudhi::Simplex_tree<>& complex)
{
	Simplex simplex;

	/// For all "ranks" of cliques (e.g. 1, 2, 3, 5, 9...)
	for (auto it = G.PH_Filtration.begin(); it != G.PH_Filtration.end(); ++it)
	{
		const int &rank = it->first;
		const auto &rankCliques = (it->second);

		/// Add cliques with size of given "rank"
		for (auto cliqueP : rankCliques)
		{
			std::vector<int> &clique = cliqueP.first;
			double weight = cliqueP.second;

			int k = (int)clique.size();

				if (k == 1)
				{
					/// Add 0-simplex
					simplex = { clique[0] };
					complex.insert_simplex_and_subfaces(simplex, 0);
				}
				else
				{
					/// Add (k-1)-simplex
					simplex = clique;
					complex.insert_simplex_and_subfaces(simplex, rank);
				}
		}
	}
}

void PH_GUDHI::buildGraphComplex(Graph& G,  Gudhi::Simplex_tree<>& complex)
{
	Simplex simplex;

	/// Add vertex:		0-simplices
	for(int v = 0; v < G.V; ++v)
	{
		simplex = { v };
		complex.insert_simplex_and_subfaces(simplex, 0);
	}

	/// Build EdgeReference
	//ReferenceArray_ini(G, true);

	/// Build WeightRank
	bool nonDecreasing = false;
	bool zeroBased = true;
	G.buildWeightRank(nonDecreasing, zeroBased);

	/// Add edges:		1-simplices
	for(int i = 0; i < G.edgesReference.size(); ++i)
	{
		int u = std::get<0>(G.edgesReference[i]);
		int v = std::get<1>(G.edgesReference[i]);
		//double w = std::get<2>(G.edgesReference[i]);
		//int edgeIndex = i;	/// "i" because edges are taken in order (nondecreasing)
		int filter = G.edgeRank(u, v);
		simplex = { u, v };
		complex.insert_simplex_and_subfaces(simplex, filter);
	}
}


void PH_GUDHI::buildCliqueComplex(Graph& G,  Gudhi::Simplex_tree<>& complex)
{
	/// IMPORT GRAPH (as 1-skeleton) FROM FILE INTO SIMPLEX-TREE COMPLEX
	/// File format: 1 simplex per line
	/// Dim1 X11 X12 ... X1d Fil1
	/// Dim2 X21 X22 ... X2d Fil2
	/// etc
	///
	/// The vertices must be labeled from 0 to n-1.
	/// Every simplex must appear exactly once.
	/// Simplices of dimension more than 1 are ignored.
	/// (*) MITXAEL: same filter can be used for two or more simplices.
	//std::string filegraph = "E:\\OneDrive\\University\\TESI\\_Resources\\Nth level - Software\\GUDHI\\GUDHI_2.1.0\\data\\filtered_simplicial_complex\\graph1_complex.fsc";
	//auto g = Gudhi::read_graph<Graph_t, Filtration_value, Vertex_handle>(filegraph);
	//complex.insert_graph(G);
	
	/// IMPORT GRAPH (as 1-skeleton) INTO SIMPLEX-TREE COMPLEX
	buildGraphComplex(G,  complex);
	
	/// EXPAND 1-skeleton UNTIL DIMENSION "max_dim" (i.e. Maximum Clique size < G.E)
	int max_dim = G.maximalCliqueSet.rbegin()->size();
	complex.expansion(max_dim);
	//std::cout << "max_dim = " << max_dim << "\n";	
}

int PH_GUDHI::intervalsFromFile(int a)
{
	///SimplexTree.- efficient and flexible data structure for representing general (filtered) simplicial complexes.
	///char* file = "E:\\OneDrive\\University\\TESI\\_Resources\\Nth level - Software\\GUDHI\\GUDHI_2.1.0\\data\\filtered_simplicial_complex\\Klein_bottle_complex.fsc";
	char* file = "E:\\OneDrive\\University\\TESI\\_Resources\\Nth level - Software\\GUDHI\\GUDHI_2.1.0\\data\\filtered_simplicial_complex\\stella_complex.fsc";
	/*	Betti numbers:	{0: 1, 1: 1}
		Betti sequence:	1 1
		Intervals & Generators:
			dim-0:
				[0.0, infinity) : [0]
			dim-1:
				[8.0, infinity) : [1,7] + [1,3] + [5,7] + [3,5]*/

	Persistence_intervals p(file);
	std::pair<double, double> min_max_ = p.get_x_range();
	std::cout << "Birth-death range : " << min_max_.first << " " << min_max_.second << std::endl;

	std::vector<double> dominant_ten_intervals_length = p.length_of_dominant_intervals(10);
	std::cout << "Length of ten dominant intervals : " << std::endl;
	for (size_t i = 0; i != dominant_ten_intervals_length.size(); ++i) {
		std::cout << dominant_ten_intervals_length[i] << std::endl;
	}

	std::vector<std::pair<double, double> > ten_dominant_intervals = p.dominant_intervals(10);
	std::cout << "Here are the dominant intervals : " << std::endl;
	for (size_t i = 0; i != ten_dominant_intervals.size(); ++i) {
		std::cout << "( " << ten_dominant_intervals[i].first << "," << ten_dominant_intervals[i].second << std::endl;
	}

	std::vector<size_t> histogram = p.histogram_of_lengths(10);
	std::cout << "Here is the histogram of barcode's length : " << std::endl;
	for (size_t i = 0; i != histogram.size(); ++i) {
		std::cout << histogram[i] << " ";
	}
	std::cout << std::endl;

	std::vector<size_t> cumulative_histogram = p.cumulative_histogram_of_lengths(10);
	std::cout << "Cumulative histogram : " << std::endl;
	for (size_t i = 0; i != cumulative_histogram.size(); ++i) {
		std::cout << cumulative_histogram[i] << " ";
	}
	std::cout << std::endl;

	std::vector<double> char_funct_diag = p.characteristic_function_of_diagram(min_max_.first, min_max_.second);
	std::cout << "Characteristic function of diagram : " << std::endl;
	for (size_t i = 0; i != char_funct_diag.size(); ++i) {
		std::cout << char_funct_diag[i] << " ";
	}
	std::cout << std::endl;

	std::vector<double> cumul_char_funct_diag =
		p.cumulative_characteristic_function_of_diagram(min_max_.first, min_max_.second);
	std::cout << "Cumulative characteristic function of diagram : " << std::endl;
	for (size_t i = 0; i != cumul_char_funct_diag.size(); ++i) {
		std::cout << cumul_char_funct_diag[i] << " ";
	}
	std::cout << std::endl;

	std::cout << "Persistence Betti numbers \n";
	std::vector<std::pair<double, size_t> > pbns = p.compute_persistent_betti_numbers();
	for (size_t i = 0; i != pbns.size(); ++i) {
		std::cout << pbns[i].first << " " << pbns[i].second << std::endl;
	}

	/// MITXAEL
	std::cout << "Boh...\n";
	///p.compute_persistent_betti_numbers();

	return 0;
}

void PH_GUDHI::barcodesIntoGraph(Graph &G, Persistent_cohomology& pcoh)
{
	/// Get intervals from barcodes
	Barcodes intervals = pcoh.get_persistent_pairs();

	for (int i = 0; i < intervals.size(); ++i)
	{
		/// Set dimension from where to get barcodes
		int dim = i;
		DataStructures::barcodeCollection barcoll;

		/// Get EndPoints from intervals, sorted
		auto endPoints = pcoh.intervals_in_dimension(dim);
		std::sort(endPoints.begin(), endPoints.end(), [](auto &left, auto &right) { return left.first < right.first; });

		/// Get Generators from barcodes
		std::vector<std::vector<int>> G_dim = {};

		if (endPoints.size() > 0)
		{
			for (auto&& pair : endPoints)
			{
				barcoll.insert(std::make_pair(std::make_pair(pair.first, pair.second), G_dim));
			}
		}

		/// Add barcodes from current dimension
		if(barcoll.size() > 0)
			G.PH_barcodes.insert(std::make_pair(dim, barcoll));
	}
}

void PH_GUDHI::holesIntoGraph(Graph &G, Persistent_cohomology& ph)
{
	
}

/*int PH_GUDHI::computePhFromSimplexTree(int a)
{

	/// TEST OF INSERTION
	std::cout << "********************************************************************" << std::endl;
	std::cout << "TEST OF INSERTION" << std::endl;
	Gudhi::Simplex_tree<> st;
	Simplex element;

	// ++ FIRST
	std::cout << "   - INSERT (0,1,2)" << std::endl;
	element = { 0, 1, 2 };
	st.insert_simplex_and_subfaces(element, 0.3);

	// ++ SECOND
	std::cout << "   - INSERT 3" << std::endl;
	element = { 3 };
	st.insert_simplex_and_subfaces(element, 0.1);

	// ++ THIRD
	std::cout << "   - INSERT (0,3)" << std::endl;
	element = { 0, 3 };
	st.insert_simplex_and_subfaces(element, 0.2);

	// ++ FOURTH
	std::cout << "   - INSERT (0,1) (already inserted)" << std::endl;
	element = { 0, 1 };
	st.insert_simplex_and_subfaces(element, 0.2);

	// ++ FIFTH
	std::cout << "   - INSERT (3,4,5)" << std::endl;
	element = { 3, 4, 5 };
	st.insert_simplex_and_subfaces(element, 0.3);

	// ++ SIXTH
	std::cout << "   - INSERT (0,1,6,7)" << std::endl;
	element = { 0, 1, 6, 7 };
	st.insert_simplex_and_subfaces(element, 0.4);

	// ++ SEVENTH
	std::cout << "   - INSERT (4,5,8,9)" << std::endl;
	element = { 4, 5, 8, 9 };
	st.insert_simplex_and_subfaces(element, 0.4);

	// ++ EIGHTH
	std::cout << "   - INSERT (9,10,11)" << std::endl;
	element = { 9, 10, 11 };
	st.insert_simplex_and_subfaces(element, 0.3);

	// ++ NINETH
	std::cout << "   - INSERT (2,10,12)" << std::endl;
	element = { 2, 10, 12 };
	st.insert_simplex_and_subfaces(element, 0.3);

	// ++ TENTH
	std::cout << "   - INSERT (11,6)" << std::endl;
	element = { 6, 11 };
	st.insert_simplex_and_subfaces(element, 0.2);

	// ++ ELEVENTH
	std::cout << "   - INSERT (13,14,15)" << std::endl;
	element = { 13, 14, 15 };
	st.insert_simplex_and_subfaces(element, 0.25);

	// Inserted simplex:
	///    1   6
	///    o---o
	///   /X\7/      4       2
	///  o---o---o---o       o
	///  2   0   3\X/8\  10 /X\
	///            o---o---o---o
	///            5   9\X/    12
	///                  o---o
	///                 11   6
	/// In other words:
	///   A facet [2,1,0]
	///   An edge [0,3]
	///   A facet [3,4,5]
	///   A cell  [0,1,6,7]
	///   A cell  [4,5,8,9]
	///   A facet [9,10,11]
	///   An edge [11,6]
	///   An edge [10,12,2]


	std::cout << "The complex contains " << st.num_simplices() << " simplices - " << st.num_vertices() << " vertices "
		<< std::endl;
	std::cout << "   - dimension " << st.dimension() << std::endl;
	std::cout << std::endl << std::endl << "Iterator on Simplices in the filtration, with [filtration value]:"
		<< std::endl;
	std::cout << "**************************************************************" << std::endl;
	std::cout << "strict graph G { " << std::endl;

	for (auto f_simplex : st.filtration_simplex_range()) {
		std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
		for (auto vertex : st.simplex_vertex_range(f_simplex)) {
			std::cout << static_cast<int>(vertex) << " -- ";
		}
		std::cout << ";" << std::endl;
	}

	std::cout << "}" << std::endl;
	std::cout << "**************************************************************" << std::endl;

	// Compute the persistence diagram of the complex
	Persistent_cohomology pcoh(st);

	/// initializes the coefficient field for homology
	int coeff_field_characteristic = 2;
	pcoh.init_coefficients(coeff_field_characteristic);
	Filtration_value min_persistence = 0.0;
	pcoh.compute_persistent_cohomology(min_persistence);

	// Output the diagram in filediag
	pcoh.output_diagram();
	return 0;
}
*/

/// Added in ".cpp" instead of ".h" because:
/// "Bottleneck.h""-> "Neighbors_finder.h" -> "CGAL/Kd_tree.h" & "CGAL/Search_traits.h" -> <mutex> = unsupported by cli/clr!!
#include "include/gudhi/Bottleneck.h"
#undef min	/// Prevent ambiguity between "windows.h" and "std::numeric_limits<>::min"
#undef max	/// Prevent ambiguity between "windows.h" and "std::numeric_limits<>::max"
double PH_GUDHI::bottleneckDistance(std::string file1, std::string file2, double tolerance = std::numeric_limits<double>::min())
{
	/// To run this program please provide as an input two files with persistence diagrams. Each file
	/// should contain a birth-death pair per line. Third, optional parameter is an error bound on the bottleneck
	/// distance (set by default to the smallest positive double value). If you set the error bound to 0, be
	/// aware this version is exact but expensive. The program will now terminate.

	/// FROM COMPLEXES
	///std::vector<std::pair<double, double>> diag1 = rips_pcoh.intervals_in_dimension(dim);
	///std::vector<std::pair<double, double>> diag2 = alpha_pcoh.intervals_in_dimension(dim);

	/// FROM FILE
	std::vector<std::pair<double, double>> diag1 = Gudhi::read_persistence_intervals_in_dimension(file1);
	std::vector<std::pair<double, double>> diag2 = Gudhi::read_persistence_intervals_in_dimension(file2);

	double b = Gudhi::persistence_diagram::bottleneck_distance(diag1, diag2, tolerance);

	return b;
}