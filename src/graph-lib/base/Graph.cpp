/// ***********************************************************
///	File	: Graph.cpp
///	About	: Graph class based on adjacency lists and many common algorithms
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
/// ******************************************************************

#include "Graph.h"				/// Graph's header file

/// Declaration of member functions

///CONSTRUCTORS

Graph::Graph()
{ 												
	this->V = 0;												/// set the given parameter
	this->adjList.resize(0);
	this->E = 0;												/// set initial value to zero
	this->K = 3;												/// set initial value to three
	this->Density = 0;											/// set initial value to density
	this->Degeneracy = 0;										/// set initial value to degeneracy
	this->Cyclomatic = 0;										/// set initial value to cyclomatic
	this->Components = 0;										/// set initial value to connected
}

Graph::Graph(int _V)
{ 											/// Constructor
	this->V = _V;												/// set the given parameter
	//this->adjList = new std::vector<std::list<std::pair<int, double>>>[_V];	/// define an array of lists
	this->adjList.resize(_V);
	this->E = 0;												/// set initial value to zero
	this->K = 3;												/// set initial value to three
	this->Density = 0;											/// set initial value to density
	this->Degeneracy = 0;										/// set initial value to degeneracy
	this->Cyclomatic = 0;										/// set initial value to cyclomatic
	this->Components = 0;										/// set initial value to connected
}

/// BASE METHODS

void Graph::reset()
{
	///this->algorithmS.clear();
	this->edgesReference.clear();
	this->cliqueSpace.clear();
	this->cycleSpace.clear();
	this->cycleSpace_w.clear();
	this->minimumCycleBasis.clear();
	this->minimumCycleBasis_w.clear();
	this->maximalCliqueSet.clear();
	this->PH_Filtration.clear();
	this->PH_barcodes.clear();
	this->vTags.clear();
	this->weightRank.clear();

	std::queue<int> empty1;	std::swap(this->BFS_spanningTree, empty1);
	std::vector<std::pair<int, int>> empty2;	std::swap(this->DFS_spanningTree, empty2);
}

bool Graph::addEdge(int src, int dst, double w)
{		
	/// Check that edge do not exist
	if(this->edgeExists(src, dst) || src == dst) return false;

	/// Push the edge's tuple in both ways (undirected ~= bidirectional)
	adjList[src].push_back(std::make_pair(dst, w));	
	adjList[dst].push_back(std::make_pair(src, w));

	///if(verbosity>2) std::cout << "Added edge: (" << src << "," << dst << "," << w << ")" << std::endl;

	E++;

	return true;
}

bool Graph::removeEdge(int src, int dst) 		// Remove an edge from Graph
{
	// first check that edge do exist
	for (auto it = adjList[src].begin(); it != adjList[src].end(); ++it)	
		// if edge exist then proceed
		if (it->first == dst) {												
			adjList[src].remove(std::make_pair(dst, (*it).second));
			adjList[dst].remove(std::make_pair(src, (*it).second));
			E--;
			return true;
		}
	//Otherwise abort deletion
	return false;
}

bool Graph::edgeExists(int src, int dst)
{
	/// For all adjacent vertices of "src"
	for (auto it = adjList[src].begin(); it != adjList[src].end(); ++it)	
	{
		/// if "dst" is adjacent to "src", then edge exists
		if (it->first == dst)												
			return true;
	}
	return false;
}

double Graph::edgeWeight(int src, int dst)
{
	/// Iterate on the vertex with less neighbours
	//if(degeneracy[src]<=degeneracy[src])
	
	if(src == dst)
	{
		return 0;
	}
	else
	{
		/// For all adjacent vertices of "src"
		for (auto it = adjList[src].begin(); it != adjList[src].end(); ++it)	
		{
			/// if "dst" is adjacent to "src", then get weight
			if (it->first == dst)
				return it->second;
		}
	}

	/// if edge do not exists
	return -1;
}

int Graph::edgeRank(int u, int v, double w)
{
	int src = (u <= v)? u : v;
	int dst = (u > v)? u : v;
	
	int rank = this->weightRank[std::make_pair(src, dst)];

	return rank;
	/*int key = 0;
	for(int i = 0; i < G.edgesReference.size(); ++i)
	{
		int t0 = std::get<0>(G.edgesReference[i]);
		int t1 = std::get<1>(G.edgesReference[i]);
		double t2 = std::get<2>(G.edgesReference[i]);

		/// If weight is -1: Look for (u,v)
		if(w == -1)
		{
			if( ( (t0==u && t1==v) || (t0==v && t1==u) ) )///&& (t2 == w) )
			{
				key = i;
				break;
			}
		}
		/// If weight is not -1: Look for (u,v,w)
		else
		{
			if( ( (t0==u && t1==v) || (t0==v && t1==u) ) && (t2 == w) )
			{
				key = i;
				break;
			}			
		}
	}
	return key;*/
}

int Graph::edgeReference(int src, int dst)
{
	for(std::size_t e = 0; e < this->edgesReference.size(); ++e)
		if( (std::get<0>(this->edgesReference[e]) == src && std::get<1>(this->edgesReference[e]) == dst)
			|| (std::get<0>(this->edgesReference[e]) == dst && std::get<1>(this->edgesReference[e]) == src) )
			return static_cast<int>(e);
	return -1;
}

int Graph::lowestWeightIndex(std::vector<int> &clique)
{
	/// Improvement: recalculate minWeight <=> removed(vertex) \in minEdge
	std::size_t k = clique.size();
	#undef max	/// Prevent ambiguity between "windows.h" and "std::numeric_limits<>::max"
	double minWgt = std::numeric_limits<double>::max();
	int minSrc, minDst;

	if(k == 1)
	{
		return 0;
	}
	else if(k == 2)
	{
		minSrc = clique.front();
		minDst = clique.back();
		minWgt = this->edgeWeight(minSrc, minDst);
		return edgeRank(minSrc, minDst);
	}
	else
	{
		for(int i = 0; i < k; ++i)
		{
			int u = clique[i];
			for(int j = i + 1; j < k; ++j)
			{
				int v = clique[j];
				double w = this->edgeWeight(u,v);
				if(this->edgeWeight(u,v) < minWgt)
				{
					minWgt = w;
					minSrc = u;
					minDst = v;
				}
			}
		}
		return edgeRank(minSrc, minDst);
	}
}

float Graph::density() {
	int V = this->V;
	int E = this->E;

	if (V > 1) {
		float e_actual = (float) E;
		float e_max = (float) (V*V-1)/2;
		return e_actual / e_max;
	}
	else
		return 0;
}
	
int Graph::degeneracy()
{
	int deg = this->E;
	for(int i = 0; i < this->V; ++i)
	{
		/// deg(v) = this->adjList[v].size();
		std::size_t deg_i = this->adjList[i].size();
		if (static_cast<int>(deg_i) < deg)
			deg = static_cast<int>(deg_i); 
	}
	return deg;
}

int Graph::cyclomaticNumber()
{
	return this->E - this->V + 1;	/// 1: n. of connected components
}

void Graph::buildGraph(boolean save2file, int n, double density) {
		
	if (verbosity > 0) std::cout << "Creating Graph based on density..." << std::endl;
	
	/// Initialize Graph
	this->V = n;
	this->E = 0;
	//this->adjList = new std::list<std::pair<int, double>>[n];
	this->adjList.resize(n);

	//printf("%.3f",density);							/// trim density to have 3 decimal places
	int e_max = ((n * (n - 1)) / 2);				/// maximum number of edges for the given number of vertices
	int e_graph = int(e_max * density);				/// expected number of edges of the Graph
	int e_vertex = (e_graph / n);					/// number of edges per vertex
	int v_incidences = (e_graph * 2 / n);			/// each edge incides "2" vertices
	int e_remainder = (e_graph % n);				/// number of remaining edges
	int dst;
	double w;

	if (e_graph >= n) {
		srand((unsigned)time(0));								/// seed random numbers to current windows time

		/// Add "V-1" consecutive edges to ensure "naive connectivity"
		for (int v = 1; v < n; ++v) {							/// Add consecutive edges to ensure Graph's connectivity
			w = (rand() % (100 + 1)) / 100.00;					/// generate random weight between 1/100 and 100/100
			int u = rand() % (v + 0);							/// generate random vertex between 0 and v
			this->addEdge(v, u, w);								/// add edge to random previous vertex
		}
			
		/// Add "e_eraph" and "e_remainder" edges
		for (int src = 0; src < n && this->E < e_graph; src++)
		{
			while ( int(this->adjList[src].size()) < v_incidences )	///while max incidences for "source" vertex not reached
			{				
				w = (rand() % (100 + 1)) / 100.00;					/// generate random weight between 1/100 and 100/100
				dst = rand() % (this->V);							/// generate random vertex between 0 and V
				if (src != dst)										/// if new edge doesn't generate self-loop
					this->addEdge(src, dst, w);						
			}
			/// If there are remainders, add ONE if possible
			if (e_remainder > 0)
			{
				bool edgeAdded = false;
				while ( (edgeAdded == false) && (static_cast<int>(this->adjList[src].size()) < v_incidences) )
				{
					int dst = rand() % (this->V);										/// generate random vertex between 0 and V
					if (src!=dst)														/// if new edge doesn't generate self-loop
						if (int(this->adjList[src].size()) < v_incidences + 1)			/// if a "remained" has not been still used for "source" vertex
							if (int(this->adjList[dst].size()) < v_incidences + 1)		/// if a "remained" has not been still used for "destination" vertex
							{
								w = (rand() % (100 + 1)) / 100.00;						/// generate random weight between 1/100 and 100/100
								edgeAdded = this->addEdge(src, dst, w);
							}	
				}
				e_remainder--;
			}
		}

		int edges = 0;
		for (int v = 0; v < this->V; v++)
			for (std::list<std::pair<int, double>>::iterator it = this->adjList[v].begin(); it != this->adjList[v].end(); ++it)
				if(v < (*it).first)
					edges++;

		this->E = edges;

		if (verbosity > 0) std::cout << "Creation completed." << std::endl;

		/// POST BUILD ACTIONS
		this->postBuild();

		/// Save graph if specified
		char buf1[4 + 1];
		snprintf(buf1, 4+1, "%04d", n);
		char buf2[5 + 1];
		snprintf(buf2, 5+1, "%.3f", density);
		std::string filename = "input_density_" + std::string(buf1) + "_" + std::string(buf2) + ".txt";
		/*std::stringstream dens;
		dens << std::setprecision(2) << density;
		std::string number1 = std::string(4 - std::to_string(V).length(), '0') + std::to_string(V);
		std::string number2 = std::string(5 - dens.str().length(), '0') + dens.str();
		std::stringstream filename;
		filename << "input_density_" << number1 << "_" << number2 << ".txt";
		//filename.str();*/
		if(save2file)	
			saveGraph(filename, 'A');
		else
			this->source = filename;
	}
	else
		std::cout << "Density is too low to ""attempt"" to generate a connected graph" << std::endl;
}

void Graph::buildGraph(boolean save2file, int n) {
	
	if (verbosity > 0) std::cout << "Creating Graph as hypercube..." << std::endl;
	
	/// initialize graph
	this->V = int(pow(2, n));							/// set number of vertices
	//this->adjList = new std::list<std::pair<int, double>>[this->V];
	this->adjList.resize(this->V);

	int e_max = ((this->V * (this->V - 1)) / 2);		/// maximum number of edges for the given number of vertices
	int e_Graph = n* static_cast<int>(pow(2, n - 1));	/// number of edges of an Qn hypercube Graph
	double density = e_max/e_Graph;						/// actual density

	/// Add "2^(n-1)" edges to every vertex
	unsigned long long max = 1ULL << n;
	for (unsigned long long vertex = 0; vertex < max; ++vertex) {
		int src = int(vertex);
		unsigned long long mask = 1;
		for (int shift_amt = 0; shift_amt < n; ++shift_amt) {
			int dst = int(vertex ^ (mask << shift_amt));
			double w = (rand() % (100 + 1)) / 100.00;	/// generate random weight between 1/100 and 100/100
			addEdge(src, dst, w);
		}
	}

	if (verbosity > 0) std::cout << "Hypercube density: " << density << std::endl;

	/// POST BUILD ACTIONS
	this->postBuild();

	/// Save graph if specified
	char buf1[4 + 1];
	snprintf(buf1, 4+1, "%04d", V);
	char buf2[5 + 1];
	snprintf(buf2, 5+1, "%05d", E);
	std::string filename = "input_hypercube_" + std::string(buf1) + "_" + std::string(buf2) + ".txt";
	/*std::string number1 = std::string(4 - std::to_string(V).length(), '0') + std::to_string(V);
	std::string number2 = std::string(4 - std::to_string(E).length(), '0') + std::to_string(E);
	std::stringstream filename;
	filename << "input_hypercube_" << number1 << "_" << number2 << ".txt";
	//filename.str()*/
	if(save2file)
		saveGraph(filename, 'A');
	else
		this->source = filename;
}

void Graph::buildGraph(boolean save2file, int n, int m)
{
	if (verbosity > 0) std::cout << "Creating Graph as Euclidean..." << std::endl;

	//if (m < n * 2.5 ) ovvero density < 0.05% ==> aided-creation
	//più grande il grafo minore la densita richiesta:
	// 4  : 4 edges
	// 5  : 7 edges
	// 7 : 8
	// 8 : 9
	// 9 : 10
	// 10 : 20 edges
	// 20 : 45 edges
	// 30 : 60 edges
	// 40 : 75
	// 50 : 80
	// 60 : 100
	// 70 : 190
	// 80 : 210
	//=> y = 10 - (x^2/6)	(for n>=10)		int q = a / b;
	//Rounded division: (a/b)+(a% b)/(b/2+b%2)
	//Rounded div: res = (a%b < (a/2.00))? a/b : (a/b)+1;
	//=> y = 12 + x/2 + (x^2)/40	(for n<10)
	// OPPURE calcolare random, e se "non-connected" => aided-creation

	/// Initialize Graph
	this->V = n;
	this->E = 0;
	//this->adjList = new std::list<std::pair<int, double>>[V];
	this->adjList.resize(V);

	/// Calculate density
	float e_actual = (float)m;
	float e_max = (float)(V*V - 1) / 2;
	float density = e_actual / e_max;

	/// Generate random coordinates for all vertices
	int max = m;
	std::vector < std::pair<double, double>> coordinates;
	srand((unsigned)time(0));										/// seed random numbers to current windows time
	for (std::size_t i = 0; i < this->V; ++i)
	{
		double x = rand() % (100 + 1) / 100.0000;					/// generate random number between 1/00 and 100/100
		double y = rand() % (100 + 1) / 100.0000;					/// generate random number between 1/00 and 100/00
		coordinates.push_back(std::make_pair(x, y));
	}

	/// Add "v-1" random but consecutive edges to ensure "naive connectivity"
	//srand((unsigned)time(0));									// seed random numbers to current windows time
	for (std::size_t i = 1; i < this->V; ++i)
	{
		int u = static_cast<int>(i);
		int v = rand() % (i);					// generate random number between 0 and i-1
		double x1 = coordinates[u].first;
		double y1 = coordinates[u].second;
		double x2 = coordinates[v].first;
		double y2 = coordinates[v].second;
		double distance = pow(x1 - x2, 2) + pow(y1 - y2, 2);
		this->addEdge(u, v, distance);
	}

	/// Add edges between random vertices, where weight=distance(u,v)
	//srand((unsigned)time(0));					// seed random numbers to current windows time
	while (this->E < max)
	{
		int u = rand() % (n);					// generate random number between 0 and n
		int v = rand() % (n);					// generate random number between 0 and n
		if (u != v)
		{
			double x1 = coordinates[u].first;
			double y1 = coordinates[u].second;
			double x2 = coordinates[v].first;
			double y2 = coordinates[v].second;
			double distance = pow(x1 - x2, 2) + pow(y1 - y2, 2);
			this->addEdge(u, v, distance);
		}
	}

	/// POST BUILD ACTIONS
	this->postBuild();

	/// Save graph if specified
	char buf1[4 + 1];
	snprintf(buf1, 4+1, "%04d", n);
	char buf2[5 + 1];
	snprintf(buf2, 5 + 1, "%05d", m);
	std::string filename = "input_random_" + std::string(buf1) + "_" + std::string(buf2) +".txt";
	/*std::string number1 = std::string(4 - std::to_string(V).length(), '0') + std::to_string(V);
	std::string number2 = std::string(4 - std::to_string(E).length(), '0') + std::to_string(E);
	std::stringstream filename;
	filename << "input_euclidean_" << number1 << "_" << number2 << ".txt";
	//filename.str()*/
	if(save2file)
		saveGraph(filename, 'A');
	else
		this->source = filename;
}

void Graph::buildGraph(std::string filename) {

	bool normalize = true;									/// shift vertices if some are missing

	/// SET INPUT FILE
	source = filename;
	filename.insert(0, ".\\input\\");
	const char *file_ = filename.c_str();
	std::ifstream fin(file_);
	if (verbosity > 0) std::cout << "Creating Graph from file " << source << std::endl;

	/// SET VARIABLES
	int V, E, src, dst;
	double w = -1;
	int MIN, MAX;
	#undef max												/// Set "max" to "std::numeric_limits<>::max"
	int lowerBound = std::numeric_limits<int>::max();
	#define max												/// Reset "max" to "windows.h"
	int upperBound = 0;

	/// READ INPUT FILE
	if (fin.is_open())
	{

		/// DETERMINE IMPORT FORMAT
		std::vector<std::string> values(3);
		std::vector<std::stringstream*> sstreams;
		std::vector<std::vector<double>> lines(3);
		int i = 0;
		while ( (values[0].empty() || values[1].empty() || values[2].empty()) && !fin.eof())
		{
			std::getline(fin, values[i]);
			sstreams.push_back(new std::stringstream(values[i]));
			i++;
		}
		fin.clear();	fin.seekg(0, std::ios::beg);			/// Return to beginning of file	
		for(i = 0; i < 3; ++i)
		{
			while (1)
			{ 
				double n;
				*sstreams[i] >> n;
				if (!*sstreams[i]) break;
				lines[i].push_back(n);
			}
		}
		char format;

		/// FORMAT A: 1st-line="V";2nd-line="e";i-th-lines="(src,dst,w)"
		if (lines[0].size() == 1 && lines[1].size() == 1 && lines[2].size() == 3)
		{
			format = 'A';

			/// ADD VERTICES
			fin >> V >> E;
			this->V = V;
			this->adjList.resize(this->V);
			
			/// ADD EDGES
			for (int i = 0; i < E; i++)
			{										
				/// ADD EDGE
				fin >> src >> dst >> w;
				if(src >= 0 && dst >= 0 && w >= 0)
				{
					this->addEdge(src, dst, w);
				}

				/// ADD VTAG
				this->vTags.insert(std::make_pair(src, std::string()));	///std::to_string(src)
				this->vTags.insert(std::make_pair(dst, std::string()));	///std::to_string(dst)

				/// TRACK BOUNDARIES
				/*MIN = min(src, dst);
				if (MIN < (lowerBound))							/// if less vertex is new then update "V"
				{
					lowerBound = MIN;
				}
				MAX = max(src, dst);
				if (MAX > (upperBound - 1))						/// if greater vertex is new then update "V"
				{
					upperBound = MAX + 1;
				}*/
			}
		}

		/// FORMAT B: i-th-lines="(src,dst,w)"
		else if (lines[0].size() == 3 && lines[1].size() == 3 && lines[2].size() == 3) 
		{
			format = 'B';

			/// COMPUTE NUMBER OF VERTICES
			while (!fin.eof())
			{									
				fin >> src >> dst >> w;

				/// TRACK BOUNDARIES
				/*MIN = min(src, dst);
				if (MIN < (lowerBound))							/// if less vertex is new then update "V"
				{
					lowerBound = MIN;
				}
				MAX = max(src, dst);
				if (MAX > (upperBound - 1))						/// if greater vertex is new then update "V"
				{
					upperBound = MAX + 1;
				}*/

				/// ADD VTAG
				this->vTags.insert(std::make_pair(src, std::string()));	///std::to_string(src)
				this->vTags.insert(std::make_pair(dst, std::string()));	///std::to_string(dst)
			}

			/// ADD VERTICES
			/*if(lowerBound != 0)									/// not zero-based
				this->V = upperBound - lowerBound;
			else												/// zero-based
				this->V = 0;*/
			MIN = this->vTags.begin()->first;		/// last on vTags
			MAX = this->vTags.rbegin()->first;		/// first on vTags
			this->V = MAX + 1;						/// index to number

			/// ADD EDGES
			this->adjList.resize(this->V);
			fin.clear();	fin.seekg(0, std::ios::beg);		/// Return to beginning of file
			while (!fin.eof())
			{											
				/// RESET VARIABLES & READ DATA
				src = dst = -1;
				w = -1.00;
				fin >> src >> dst >> w;

				/// ADD EDGES (IF VALID)
				if(src >= 0 && dst >= 0 && w >= 0)
				{
					this->addEdge(src, dst, w);									
				}
			}
		}

		/// FORMAT C: 1st-line="v"; i-th-lines="(src,w-1,w-2,...w-n)"
		else if(lines[0].size() == 1 && lines[1].size() > 1 && lines[2].size() > 1)
		{
			format = 'C';
		
			/// SET VERTICES
			fin >> V;
			this->V = V;

			/// ADD EDGES
			///TODO: Should add to a matrix (std::vector<std::vector<int>>)
			this->adjList.resize(this->V);
			fin.clear();	fin.seekg(0, std::ios::beg);		/// Return to beginning of file
			while (!fin.eof())
			{											
				/// RESET VARIABLES & READ DATA
				src = dst = -1;
				w = -1.00;
				fin >> src >> dst >> w;

				/// ADD EDGES (IF VALID)
				if(src >= 0 && dst >= 0 && w >= 0)
				{
					this->addEdge(src, dst, w);

					/// ADD VTAG
					this->vTags.insert(std::make_pair(src, std::string()));	///std::to_string(src)
					this->vTags.insert(std::make_pair(dst, std::string()));	///std::to_string(dst)
				}
			}
		}

		/// FORMAT UNKNOWN
		else  /// if(values.size() > 3)
		{
			format = 'U';
			if (verbosity > 0) std::cout << "Unknown file format. Unable to create graph." << std::endl;
		}

		/// CLOSE INPUT FILE
		fin.close();

		/// NORMALIZE VERTICES IF NEEDED
		if(normalize == true && this->V != (int)this->vTags.size())
		{
			int newVertex = 0;
			int offset = 0;
			DataStructures::verticesTags::iterator vtag_it = this->vTags.begin();
			///for (auto it = this->vTags.cbegin(); it != this->vTags.cend(); /* no increment to handle deletions*/)
			while( newVertex < (int)this->adjList.size() )
			{
				int oldVertex = (*vtag_it).first;

				/// IF NEW IS DIFFERENT FROM OLD THEN SHIFT DOWN (i.e. REMOVE NEW VERTEX FROM ADJ LIST)
				if( (newVertex + offset) < oldVertex )
				{
					/// REMOVE OLD VERTEX (SHIFT DOWN)
					this->adjList.erase( this->adjList.begin() + newVertex );				

					/// INCREASE "offset" BECAUSE OF NEWVERTEX
					++offset;
				}

				/// IF NEW MATCHES OLD THEN UPDATE
				else
				{
					/// SAVE OLD VALUE AS NEW TAG
					std::swap(this->vTags[newVertex], std::to_string(oldVertex));

					/// REMOVE OLD TAG & MOVE TO NEXT "vTag"
					this->vTags.erase(vtag_it++);	

					/// REPLACE OLD WITH NEW ON ADJACENCY LIST
					for (auto& list : this->adjList)
					{
						for (auto& element : list)
						{
							if(element.first == oldVertex)
								element.first = newVertex;
						}
					}

					/// INCREASE "newVertex"
					++newVertex;
				}
			}

			/// UPDATE GRAPH.VERTICES
			this->V = (int) this->adjList.size();
		}

		/// POST BUILD ACTIONS
		this->postBuild();
	}
	else
	{
		if (verbosity > 0) std::cout << "The file either cannot be read or doesn't exist." << std::endl;
	}
}
	
void Graph::buildGraph()
{	
	std::cout << "Select a Graph from menu..." << std::endl;

	/// PREPARE
	const int size = sizeof(WIN32_FIND_DATA);
	int option = 0;
	const char *source = ".\\input\\*.txt";
	WCHAR target[size];
	MultiByteToWideChar(CP_ACP, 0, source, -1, target, size);
	LPCWSTR file_ = target;

	/// SHOW FILES IN MENU
	std::vector<std::wstring> files;	//container of files (related to an index)
	WIN32_FIND_DATA search_data;
	memset(&search_data, 0, sizeof(WIN32_FIND_DATA));
	HANDLE handle = FindFirstFile(file_, &search_data);
	do
	{
		option++;
		files.push_back(search_data.cFileName);
		std::wcout << option << ": " << (search_data.cFileName) << std::endl;	//printf("Found file: %ws\r\n", search_data.cFileName);
	} while ((handle != INVALID_HANDLE_VALUE) && (FindNextFile(handle, &search_data) != FALSE));

	/// SELECT OPTION FROM MENU
	option = 0;
	do {
		std::cin >> option;
		if ((option > 0) && (option <= files.size())) {
			std::string fname(files[option - 1].begin(), files[option - 1].end());
			buildGraph(fname);
		}
		else
			std::cout << "Please select a valid option...";
	} while ((option <= 0) || (option > files.size()));

	/// POST BUILD ACTIONS
	this->postBuild();
}

void Graph::saveGraph(std::string filename, char format)
{	
	if (verbosity > 0) std::cout << "Saving Graph to file..." << std::endl;
	
	std::string fn_original = filename;
	filename.insert(0, ".\\input\\");
	const char *file_ = filename.c_str();

	std::ofstream fout(file_);
	if (fout.is_open()) {
		
		/// Format=Adjacency
		if (format == 'A')
		{
			fout << this->V << "\n";
			fout << this->E << "\n";
		}

		/// Format=Matrix-Market
		if (format == 'M')
		{
			fout << "%%MatrixMarket matrix coordinate real symmetric" << "\n";
			fout << "%" << "\n";
			fout << "%  matrix: mxn matrix" << "\n";
			fout << "%  coordinate : only(coordinates of) nonzero entries are provided." << "\n";
			fout << "%  real : arithmetic field of values is Real" << "\n";
			fout << "%  symmetric : diagonal entries are zero(and therefore omitted).Only entries" << "\n";
			fout << "%  in the lower triangular portion need be supplied." << "\n";
			fout << "%  general : matrix format is not taking advantage of any symmetry properties." << "\n";
			fout << "%" << "\n";
			fout << "%  1st - line : rows cols edges" << "\n";
			fout << "%  2nd - onwards : src dst weight" << "\n";
			fout << "%  Conditions : vertices start at 1, no loops." << "\n";
			fout << "%" << "\n";
			fout << this->V << " " << this->V << " " << this->E << "\n";
		}

		for (auto v = 0; v < this->V; ++v)
		{
			for (auto e = adjList[v].begin(); e != adjList[v].end(); ++e)
			{
				/// prevent repeating edges (a->b and b->a)
				if (v < (*e).first)					
				{
					int src = v;
					int dst = (*e).first;
					double weight = (*e).second;	// *100.00;
					fout << src << " " << dst << " " << weight << "\n";
				}
			}
		}
		fout.close();

		/// Keep source filename in graph
		this->source = fn_original;
	}
}

void Graph::postBuild()
{
	this->Density = this->density();
	this->Degeneracy = this->degeneracy();
	this->Cyclomatic = this->cyclomaticNumber();
	this->Components = this->isConnected("bfs");
}

bool Graph::isConnected(std::string mode) {

	if (mode == "bfs") {
		this->bfs(0);
		if (BFS_spanningTree.size() == this->V)	/// all vertices must be spanned in the tree to state that G is connected
			return true;
		else
			return false;
	}

	if (mode == "dfs") {
		int *color = new int[V];
		int *parents = new int[V];
		this->dfs(0, 0, color, parents);
		for (std::size_t i = 0; i < this->V; ++i)	/// all vertices must be visited to state that G is connected
			if (color[i] == 0)
				return false;
			else
				return true;
	}

	return false;
}

/// GRAPH ALGORITHMS

void Graph::bfs(int root)
{
	/// INITIALIZE
	std::queue<int> T;				/// Spanning Tree
	std::queue<int> q;				/// Queue for BFS
	int *color;						/// 0=white: not visited; 1=gray: visited; 2=black: completed
	color = new int[V];
	for (int i = 0; i < V; i++)
		color[i] = 0;
	
	/// Choose a starting vertex
	int src = root;
	//srand((unsigned)time(0));							/// seed random numbers to current windows time
	//src = rand() % (V + 0);							/// generate random vertex between 0 and V-1

	/// Clear BFS
	std::queue<int> empty;
	std::swap(BFS_spanningTree, empty);

	/// Mark vertex
	color[src] = 1;
	/// Add vertex to queue
	q.push(src);
	/// Add vertex to Tree
	T.push(src);

	/// COMPUTE SpanningTtree
	while (!q.empty())
	{
		/// Choose the vertex in fron of the queue
		int u = q.front();
		q.pop();	

		std::list<std::pair<int, double>>::iterator it;
		for (it = adjList[u].begin(); it != adjList[u].end(); ++it)	/// Recur for all vertices adjacent to this vertex
		{
			if (color[(*it).first] == 0)
			{
				/// Mark vertex
				color[(*it).first] = 1;
				/// Add vertex to queue
				q.push((*it).first);
				/// Add vertex to Tree
				T.push((*it).first);
			}
		}
		color[u] = 2;
	} 

	BFS_spanningTree = T;

	/*while (!T.empty())	//print T
	{
		int item = T.front();
		std::cout << item  << ",";
		T.pop();
	}*/
}

bool Graph::dfs(int v, int parent, int color[], int parents[])
{
	// Mark the current node as visited
	color[v] = 1;
	parents[v] = parent;

	// Recur for all the vertices adjacent to this vertex
	std::list<std::pair<int, double>>::iterator it;
	for (it = adjList[v].begin(); it != adjList[v].end(); ++it)
	{
		if (color[(*it).first] == 0)					// If an adjacent is not visited, then recur "dfs" for it
		{
			if (dfs((*it).first, v, color, parents))
				return true;
		}
		else if ((*it).first != parent)				// Else, if an adjacent is visited and not parent of current vertex, then there is a cycle.
		{
			return true;
		}
	}
	color[v] = 2;								// Mark the current node as completed

	return false;
}

void Graph::bfsFast()
{
	/// INITIALIZE
	std::queue<int> q;				// Queue for BFS
	int *visited = new int[V];		// 0=white: not visited; 1=gray: visited; 2=black: completed
	for (int i = 0; i < V; i++)
		visited[i] = 0;
	srand((unsigned)time(0));							// seed random numbers to current windows time
	int src = rand() % (V + 0);							// generate random weight between 0 and V-1

	/// COMPUTE SpanningTtree
	q.push(src);
	do
	{
		int u = q.front();
		q.pop();

		std::list<std::pair<int, double>>::iterator it;
		for (it = adjList[u].begin(); it != adjList[u].end(); ++it)	// Recur for all vertices adjacent to this vertex
		{
			if (visited[(*it).first] == 0)
			{
				q.push((*it).first);
				visited[(*it).first] = 1;
				(u<(*it).first)? DFS_spanningTree.push_back(std::make_pair(u, (*it).first)) : DFS_spanningTree.push_back(std::make_pair((*it).first,u));
			}
		}
		visited[u] = 2;
	} while (!q.empty());
}

void Graph::dfsFast_caller()
{
	srand((unsigned)time(0));								// seed random numbers to current windows time
	int v = rand() % (V + 0);								// generate random weight between 0 and V-1
	int *visited = new int[V];								// define array for vertices check
	for (int i = 0; i < V; i++) visited[i] = 0;				// initialize to zero
	dfsFast(v, visited);									// start dfsFast
}

void Graph::dfsFast(int v, int visited[])
{
	// Mark the current node as visited
	visited[v] = 1;

	// Recur for all the vertices adjacent to this vertex
	std::list<std::pair<int, double>>::iterator it;
	for (it = adjList[v].begin(); it != adjList[v].end(); ++it)
	{
		if (visited[(*it).first] == 0)
		{
			dfsFast((*it).first, visited);
			(v<(*it).first) ? DFS_spanningTree.push_back(std::make_pair(v, (*it).first)) : DFS_spanningTree.push_back(std::make_pair((*it).first, v));
		}
	}
}

/// CYCLES CONTROL

bool Graph::existCycle(void)
{
	int *color = new int[V];						// Mark all the vertices as not visited and not part of recursion stack
	int *parents = new int[V];
	for (int i = 0; i < V; i++)
		color[i] = 0;

	for (int u = 0; u < V; u++)						// Call the recursive helper function to detect cycle in different DFS trees
	{
		if (color[u] == 0)							// Recur for "u" only if it has not been visited
			if (dfs(u, -1, color, parents))
			{
				std::cout << "Graph contains cycle(s)\n";

				for (int j = 0; j < V; j++)
					if (color[j] == 1)
					std::cout << j << " ";
				std::cout << std::endl;

				return true;
			}
	}
	std::cout << "Graph doesn't contain cycle\n";

	return false;
}

/*int Graph::addCycleEdge(int src, int dst) {
	std::size_t lastPos = this->cycleSpace.size() - 1;							// position of the last entered cycle
	for (std::size_t i = 0; i < this->cycleSpace[lastPos].size(); ++i)			// travel all edges, enabling them if matching the edge
		if (((this->edgesReference[i].first == src) && (this->edgesReference[i].second == dst)) ||
			((this->edgesReference[i].first == dst) && (this->edgesReference[i].second == src))) {
			this->cycleSpace[lastPos][i] = 1;
			return static_cast<int>(i);
			/// Add weight now! to optimize the "cycleWeight" and "sorting" functions
			//for (std::list<std::pair<int, double>>::iterator it = G1.adjList[src].begin(); it != G1.adjList[src].end(); it++)
			//if ((*it).first == dst)
			//G1.cycleSpace_w[lastPos] += (*it).second;
		}
	return 0;
}*/

void Graph::addCycleToCS(std::vector<ut> &cycle, double weight)
{
	this->cycleSpace.push_back(cycle);

	if(weight >= 0)
		this->cycleSpace_w.push_back(weight);
	else
		this->cycleSpace_w.push_back(this->cycleWeight(static_cast<int>(this->cycleSpace.size())-1, cycle));
}

void Graph::removeCycleFromCS(int pos)
{
	if(this->cycleSpace.size() > pos)
	{
		this->cycleSpace.erase(this->cycleSpace.begin() + pos);
		this->cycleSpace_w.erase(this->cycleSpace_w.begin() + pos);
	}
}

void Graph::addCycleToMCB(int pos)
{
	if(this->cycleSpace.size() > pos)
	{
		this->minimumCycleBasis.push_back(this->cycleSpace[pos]);
		this->minimumCycleBasis_w.push_back(this->cycleSpace_w[pos]);	
	}
}

void Graph::addHoleToMCB(std::vector<ut> &hole, double weight)
{
	this->minimumCycleBasis.push_back(hole);

	if(weight >= 0)
		this->minimumCycleBasis_w.push_back(weight);
	else
		this->minimumCycleBasis_w.push_back(this->cycleWeight(static_cast<int>(this->minimumCycleBasis.size())-1, hole));
}

void Graph::removeCycleFromMCB(int pos)
{
	if(this->minimumCycleBasis.size() > pos)
		this->minimumCycleBasis.erase(this->minimumCycleBasis.begin() + pos);
	if(this->minimumCycleBasis_w.size() > pos)
		this->minimumCycleBasis_w.erase(this->minimumCycleBasis_w.begin() + pos);
}

double Graph::cycleWeight(int pos, std::vector<ut> &cycle)
{
	int src, dst;
	double cw = 0;

	if(!this->cycleSpace_w.empty())
	{
		if ( cycle.size() == this->cycleSpace_w.size() )
			cw = this->cycleSpace_w[pos];
		if ( cycle.size() == this->minimumCycleBasis_w.size() )
			cw = this->minimumCycleBasis_w[pos];
	}
	else
	{
		for (std::size_t i = 0; i < cycle.size(); i++) {
			src = std::get<0>(this->edgesReference[i]);
			dst = std::get<1>(this->edgesReference[i]);
			if (cycle[i] == 1)
			{
				for (std::list<std::pair<int, double>>::iterator it = this->adjList[src].begin(); it != this->adjList[src].end(); ++it) {
					if ((*it).first == dst)
						cw = cw + (*it).second;
				}
			}
		}
	}

	return cw;
}

double Graph::minimumCycleBasisWeight() {

	double minimumCycleBasis_single_weight;
	double minimumCycleBasis_total_weight = 0;
	for (std::size_t r = 0; r < this->minimumCycleBasis.size(); ++r) {
		minimumCycleBasis_single_weight = cycleWeight(static_cast<int>(r), this->minimumCycleBasis[r]);
		minimumCycleBasis_total_weight += minimumCycleBasis_single_weight;
	}
	if (verifyMCB() == true)
		return minimumCycleBasis_total_weight;
	else
		return ((minimumCycleBasis_total_weight) * (-1));
}

bool Graph::verifyMCB() {

	if (this->minimumCycleBasis.size() > 0) {
			
		/// Create a Result vector of size "m" (m: incident edges of G)
		std::vector<ut> result(this->minimumCycleBasis[0].size(), 0);
			
		/// Combine all cycles (rows) in MCB using the "symmetric difference"
		for (std::size_t i = 0; i < minimumCycleBasis.size(); ++i)
			result = symmetricDifference_incidence(result, this->minimumCycleBasis[i]);

		/// Check the consistency of the Result vector
		int consistent = result[0];

		/// Look for at least one 0 in the Result vector
		for (std::size_t i = 0; i < result.size(); ++i)
			consistent = consistent & result[i];  // AND: Always 1 unless there's at least one 0 
		/// If the result vector is consistent (i.e. contains all 1's) then all edges have been covered
		
		if (verbosity > 2) std::cout << std::endl << "MBC consistency = " << consistent << std::endl;
		return (consistent == 1)? true : false;
	}
	else
		return false;
}

/// DATA VISUALIZATION

void Graph::displayGraph() {
	std::cout << std::endl;
	for (int v = 0; v < this->V; v++)
	{
		std::cout << "" << v << " : ";
		std::list<std::pair<int, double>>::iterator it;
		for (it = this->adjList[v].begin(); it != this->adjList[v].end(); ++it)
		{
			std::cout << (*it).first << "(" << (*it).second << ") ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void Graph::displayCycleSpace() {

	std::cout << "The space of candidate cycles has " << this->cycleSpace.size() << " elements:" << std::endl;

	/// Display the binary matrix
	if (verbosity > 2)
	{
		for (std::size_t e = 0; e < this->edgesReference.size(); ++e)
			std::cout << "(" << std::get<0>(this->edgesReference[e]) << "," << std::get<1>(this->edgesReference[e]) << ")";
		std::cout << std::endl;

		for (std::size_t r = 0; r < this->cycleSpace.size(); ++r) {
			for (std::size_t c = 0; c < this->cycleSpace[r].size(); ++c)
				std::cout << (int)this->cycleSpace[r][c] << " ";
			std::cout << " ... " << this->cycleSpace_w[static_cast<int>(r)] << std::endl;
		}
		std::cout << std::endl;
	}

	/// Display the cicles and edges
	if (verbosity > 1)
	{
		for (std::size_t r = 0; r < this->cycleSpace.size(); ++r) {
			std::cout << "C" << r << " = { ";
			for (std::size_t c = 0; c < this->cycleSpace[r].size(); ++c) {
				if (this->cycleSpace[r][c] == 1)
					std::cout << "(" << std::get<0>(this->edgesReference[c]) << "," << std::get<1>(this->edgesReference[c]) << ") ";
			}
			double minimumCycleBasis_single_weight = this->cycleSpace_w[static_cast<int>(r)];
			std::cout << " }  (w = " << minimumCycleBasis_single_weight << ")" << std::endl;
		}
		std::cout << std::endl;
	}
}

void Graph::displayPath(std::vector<ut> pathReference)
{
	/// Display the full path
	std::cout << "Path = { ";
	for (std::size_t i = 0; i < pathReference.size(); ++i) {
		if (pathReference[i] == 1)
			std::cout << "(" << std::get<0>(this->edgesReference[i]) << "," << std::get<1>(this->edgesReference[i]) << ") ";
	}
	std::cout << std::endl;
}

double Graph::displayMCB(bool verify) {
		
	int minimumCycleBasisSize = int(this->minimumCycleBasis.size());
	double minimumCycleBasis_total_weight = 0;

	/// If a minimum cycle basis (MCB) was found
	if (verify==false || (minimumCycleBasisSize == (this->E - this->V + 1)) )
	{
		//SORT MCB
		/*for (std::size_t gap = this->minimumCycleBasis.size() / 2; gap > 0; gap /= 2)
		{
			for (int i = gap; i < this->minimumCycleBasis.size(); ++i)
				for (int j = i - gap; j >= 0 && cycleWeight(j, this->minimumCycleBasis) > cycleWeight(j + gap, this->minimumCycleBasis); j -= gap) {
					std::swap(this->minimumCycleBasis[j], this->minimumCycleBasis[j + gap]);
				}
		}*/

		/// Display the binary matrix
		if (verbosity > 0)
		{
			std::cout << "The binary matrix of the minimum cycle basis (MCB) is: " << std::endl;
			for (std::size_t r = 0; r < minimumCycleBasis.size(); ++r) {
				for (std::size_t c = 0; c < minimumCycleBasis[r].size(); ++c)
					std::cout << (int)minimumCycleBasis[r][c] << " ";
				std::cout << std::endl;
			}

			//Display the cicles and edges
			std::cout << std::endl << "The minimum cycle basis (MCB) is: " << std::endl;
			for (std::size_t r = 0; r < minimumCycleBasisSize; ++r) {
				//std::cout << "{ ";
				std::cout << "C" << r << " = { ";
				for (std::size_t c = 0; c < this->minimumCycleBasis[r].size(); ++c) {
					if (this->minimumCycleBasis[r][c] == 1)
						std::cout << "(" << std::get<0>(this->edgesReference[c]) << "," << std::get<1>(this->edgesReference[c]) << ") ";
				}
				double minimumCycleBasis_single_weight = this->minimumCycleBasis_w[static_cast<int>(r)];
				minimumCycleBasis_total_weight += minimumCycleBasis_single_weight;
				std::cout << " }  (w = " << minimumCycleBasis_single_weight << ")" << std::endl;
			}
			if (verbosity > 0) std::cout << "The total weight of the MCB is: " << minimumCycleBasis_total_weight << std::endl;
		}
		else
		{
			for (std::size_t r = 0; r < minimumCycleBasisSize; ++r)
				minimumCycleBasis_total_weight += this->minimumCycleBasis_w[static_cast<int>(r)];
			std::cout << "The total weight of the MCB is: " << minimumCycleBasis_total_weight << std::endl;
		}
		
	}

	/// Otherwise, no minimum cycle basis (MCB) was found
	return minimumCycleBasis_total_weight;
}

int Graph::displayCliques(bool sorted, int Kmin)
{
	int cliqueN = 0;

	/// Sort cliques from smallest to largest
	std::sort(this->maximalCliqueSet.begin(),this->maximalCliqueSet.end(), [](const std::vector<int> & a, const std::vector<int> & b){ return a.size() < b.size(); });

	/// Display clique
	
	std::cout << "Maximal Cliques (k>=" << Kmin << "): " << this->maximalCliqueSet.size() << std::endl;

	for (auto clique : this->maximalCliqueSet)
	{
		int cliqueSize = static_cast<int>(clique.size());
		if(cliqueSize >= Kmin)
		{
			cliqueN++;
			if (verbosity > 0)
			{
				std::cout << "Clique of size " << cliqueSize << ": [";
				for (int node : clique)
				{
					std::cout << node << " ";
				}
				std::cout << "]" << std::endl;
			}
		}
	}
	return cliqueN;
}

/// DATA ORDERING

void Graph::shellSortCycles(int n) {

	//static constexpr int gap_sequence[] = { 13, 9, 5, 2, 1 };
	//static constexpr int n = G.cycleSpace.size() / G.cycleSpace[0].size();

	for (int gap = n / 2; gap > 0; gap /= 2)
		//for (int gap : gap_sequence) if (gap < n)
	{
		for (int i = gap; i < n; ++i)
			for (int j = i - gap; j >= 0 && this->cycleSpace_w[j] > this->cycleSpace_w[j + gap]; j -= gap) {
				std::swap(this->cycleSpace[j], this->cycleSpace[j + gap]);
				std::swap(this->cycleSpace_w[j], this->cycleSpace_w[j + gap]);
				/*std::vector<int> tmp = G.cycleSpace[j];
				G.cycleSpace[j] = G.cycleSpace[j+gap];
				G.cycleSpace[j+gap] = tmp;*/
				/*for (std::size_t i = 0; i < G.cycleSpace[p].size(); i++) {
				int tmp = G.cycleSpace[j][i];
				G.cycleSpace[j][i] = G.cycleSpace[j+gap][i];
				G.cycleSpace[j+gap][i] = tmp;
				}*/
			}
	}
}

void Graph::quickSortCycles(int p, int r) {

	if (p < r) {					//(G, 0, cyclespace.size())
		int j = partition(p, r);
		quickSortCycles(p, j - 1);
		quickSortCycles(j + 1, r);
	}

}

int Graph::partition(int p, int r) {
	double pivot = this->cycleSpace_w[r];
	while (p < r)
	{
		while (this->cycleSpace_w[p] < pivot)
			p++;

		while (this->cycleSpace_w[r] > pivot)
			r--;
		if (this->cycleSpace_w[p] == this->cycleSpace_w[r])
			p++;
		else if (p < r)
		{
			std::swap(this->cycleSpace[p], this->cycleSpace[r]);
			std::swap(this->cycleSpace_w[p], this->cycleSpace_w[r]);
			/*for (std::size_t i = 0; i < G.cycleSpace[p].size(); i++) {
			int tmp = G.cycleSpace[p][i];
			G.cycleSpace[p][i] = G.cycleSpace[r][i];
			G.cycleSpace[r][i] = tmp;
			}*/
		}
	}

	return r;
}

template <typename T>
void Graph::quickSortPlus(T * begin, T * end)
{
	std::size_t len = std::distance(begin, end);
	if (len <= 1) return;
	T * pivot = end - 1;
	T * mid = std::partition(begin, pivot,
	[=](T value){ return value < *pivot; });
	std::swap(*mid, *pivot);
	std::sort(begin, mid);
	std::sort(mid + 1, end);

	///Usage:
	/// double* ptrToFirst = &this->cycleSpace_w.front();
	/// double* ptrToLast = &this->cycleSpace_w.back();
	/// this->quickSortPlus(ptrToFirst, ptrToLast);
}

void Graph::defaultSort(std::vector<std::vector<ut>>& cycle_edges, std::vector<double>& cycle_weights)
{
	//std::sort(cycle_edges.begin(), cycle_edges.end(), [&cycle_weights](auto a, auto b) {return a < b;});
	//std::sort(cycle_weights.begin(), cycle_weights.end(), [](auto a, auto b) {return a < b;});

	/// Create a vector of "common" indexes
	std::vector<std::size_t> index_vec;
	for (std::size_t i = 0; i != cycle_edges.size(); ++i)
		index_vec.push_back(i);

	/// Sort the vector of "common" indexes
	std::sort(index_vec.begin(), index_vec.end(),
		[&](std::size_t a, std::size_t b)
			{ return cycle_weights[a] < cycle_weights[b]; }
	);

	/// In-place update based on "common" indexes
	double temp_ce;
	std::vector<ut> temp_cw;
	for(size_t i = 0; i < cycle_weights.size(); ++i)
	{
		size_t j, k;
		if(i != index_vec[i])
		{
			temp_ce = cycle_weights[i];
			temp_cw = cycle_edges[i];
			k = i;
			while(i != (j = index_vec[k]))
			{
				cycle_weights[k] = cycle_weights[j];
				cycle_edges[k] = cycle_edges[j];
				index_vec[k] = k;
				k = j;
			}
			cycle_weights[k] = temp_ce;
			cycle_edges[k] = temp_cw;
			index_vec[k] = k;
		}
	}
}

void Graph::sortCycleSpace() {
	int size_CS = (this->E - this->V + 1);

	/// if the array is fairly small
	if (size_CS < 250)
	{
		this->quickSortCycles(0, static_cast<int>(this->cycleSpace.size())-1);
	}
	/// if the array is large
	else				
	{
		this->shellSortCycles(static_cast<int>(this->cycleSpace.size()));
	}

	if (verbosity > 1) std::cout << "    Cycle space successfully sorted (by weight) in non-decreasing order." << std::endl;
}

/// ANCILLARY TOOLS

void Graph::lexicalAdjList()
{
	for(std::size_t i = 0; i < adjList.size(); ++i)
	{
		std::vector<std::pair<int, double>> temp; 
		for (auto it = adjList[i].begin(); it != adjList[i].end(); ++it)
		{
			temp.push_back(*it);
			adjList[i].pop_front();
		}

		std::sort(temp.begin(), temp.end());

		/*for (auto it = temp.begin(); it != temp.end(); ++it)
		{
			adjList[i].insert(*it);
		}*/
	}
}

void Graph::buildWeightRank(bool nondecreasing, bool zeroBased)
{
	/// 0: reserve "0" to vertices and start from "1".	-1: start from "0" as vertices.
	int idx = (zeroBased)? -1 : 0;
	///TODO: FindCliques() has index problems if 1-based

	///Ascending order
	if(nondecreasing)
	{
		#undef min															/// Prevent ambiguity between "windows.h" and "std::numeric_limits<>::min"
		const double min_double = -1;//std::numeric_limits<double>::epsilon();
		double lastw = min_double;

		for (std::vector<std::tuple<int, int, double>>::iterator it = this->edgesReference.begin(); it != this->edgesReference.end(); ++it)
		{
			int v1 = std::get<0>(*it);
			int v2 = std::get<1>(*it);
			double ew = std::get<2>(*it);

			if(ew > lastw) {
				idx++;
				lastw = ew;
			}

			std::pair<int,int> edge = std::make_pair(v1, v2);
			this->weightRank.insert(std::make_pair(edge, idx));
		}
	}

	/// Descending order
	else
	{
		#undef max															/// Prevent ambiguity between "windows.h" and "std::numeric_limits<>::max"
		const double max_double = std::numeric_limits<double>::max();
		double lastw = max_double;
		//for (int i = this->edgesReference.size(); i --> 0; )
		for (std::vector<std::tuple<int, int, double>>::reverse_iterator it = this->edgesReference.rbegin(); it != this->edgesReference.rend(); ++it)
		{
			int v1 = std::get<0>(*it);
			int v2 = std::get<1>(*it);
			double ew = std::get<2>(*it);

			if(ew < lastw) {
				idx++;
				lastw = ew;
			}

			std::pair<int,int> edge = std::make_pair(v1, v2);
			this->weightRank.insert(std::make_pair(edge, idx));
		}
	}
}

bool ComparePG::operator()(const std::pair<int, double>& a, const std::pair<int, double>& b) {
		return  a.second>b.second; //return true if 1st pair's weight is GREATER than 2nd pair's weight
}

bool ComparePL::operator()(const std::pair<int, double>& a, const std::pair<int, double>& b) {
	return  a.second<b.second; ///return true if 1st pair's weight is LESS than 2nd pair's weight
}

bool ComparePFL::operator()(const std::pair<int, double>& a, const std::pair<int, double>& b) {
	return  a.first<b.first; ///return true if 1st pair's index is LESS than 2nd pair's index
}

bool CompareTL::operator()(const std::tuple<int, int, double> &lhs, const std::tuple<int, int, double> &rhs) {
		return std::get<2>(lhs) < std::get<2>(rhs); //return true if 1st pair's weight is LESS than 2nd tuple's weight
}

bool CompareTR::operator()(const std::tuple<int, int, double> &lhs, const std::tuple<int, int, double> &rhs) {
		return std::get<2>(lhs) > std::get<2>(rhs); //return true if 1st pair's weight is GREATER than 2nd tuple's weight
}

size_t HashPair::operator()(const std::pair<std::vector<int>, double>& p)
{
	int vhash;
	std::hash<int> hash_i;
	for (const int& elem : p.first)
		vhash ^= hash_i(elem);
	return vhash ^ (unsigned&)p.second;
	//return p.first.size() ^ (unsigned&)p.second;
}

std::vector<ut> Graph::symmetricDifference_incidence(std::vector<ut> V1, std::vector<ut> V2) {
	std::vector<ut> result(V1.size());
	if (V1.size() == V2.size())
	{
		for (std::size_t i = 0; i < result.size(); ++i)
			//result[i] = V1[i] ^ V2[i];		/// V1 XOR V2 : 1 iif !=
			result[i] = (V1[i] + V2[i] == 0)? 0 : 1;
	}
	if (verbosity > 2)
	{
		std::cout << std::endl << "Symmetric difference: ";
		for (std::size_t j = 0; j < result.size(); ++j)
			std::cout << result[j] << " ";
	}
	return result;
}

std::vector<std::tuple<int, int, double>> Graph::symmetricDifference_edges(std::vector<std::tuple<int, int, double>> V1, std::vector<std::tuple<int, int, double>> V2) {
	std::vector<std::tuple<int, int, double>> a = V1;	// { std::make_pair(1,3), std::make_pair(2,3), std::make_pair(4,5) };
	std::vector<std::tuple<int, int, double>> b = V2;	// { std::make_pair(1,3), std::make_pair(3,4), std::make_pair(4,5) };	
		
	std::vector<std::tuple<int, int, double>> SD(V);														// 0  0  0  0  0  0  0  0  0  0
	std::vector<std::tuple<int, int, double>>::iterator it;

	std::sort(a.begin(), a.end(), CompareTL());													//  5 10 15 20 25
	std::sort(b.begin(), b.end(), CompareTL());													// 10 20 30 40 50

	it = std::set_symmetric_difference(a.begin(), a.end(), b.begin(), b.end(), SD.begin());		//  5 15 25 30 40 50  0  0  0  0
	SD.resize(it - SD.begin());																	//  5 15 25 30 40 50

	if (verbosity > 2) {
		for (it = SD.begin(); it != SD.end(); ++it)
			std::cout << "(" << std::get<0>((*it)) << "," << std::get<1>((*it)) << ")";
		std::cout << '\n';
	}

	return SD;
}

/*std::vector<int> symmetricDifference(Graph &G1, int a, int b) {
	std::vector<int> SD(G1.cycleSpace[a].size());

	for (std::size_t i = 0; i < G1.cycleSpace[a].size(); i++)
		if (G1.cycleSpace[a][i] != G1.cycleSpace[b][i])
			SD[i] = 1;
		else
			SD[i] = 0;
	return SD;
}*/

/*bool Graph::AncestorCheck(int u, int v) //Check if "u" is the ancestor of "v"
{
	if (previsit(x) < previsit(y) && postvisit(x) > postvisit(y))
		return 1;
	return 0;
}*/

/*int Graph::LCACheck(int u, int v) //Check if "u" is the ancestor of "v"
{
	int *color = new int[V];
	int *parents = new int[V];
	dfs(u, -1, color, parents);

	if (AncestorCheck(u, v))
		then u is the LCA of v
	return;
	if (AncestorCheck(v, u))
		then v is the LCA of u
	return;
	Let a node temp = v;
	while (!(AncestorCheck(temp, u) && AncestorCheck(temp, v))
		temp = parent of temp;
	return temp;
}*/