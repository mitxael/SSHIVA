/// ******************************************************************
///	File	: Javaplex.cpp
///	About	: Javaplex PH methods
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// -  https://docs.oracle.com/javase/7/docs/technotes/guides/jni/spec/functions.html
/// ******************************************************************

#include "Javaplex.h"

JVM JVMinstance;

HINSTANCE hinstLib;
HMODULE jvmHandle;
JavaVM *jvm(0);								/// Pointer to the JVM (Java Virtual Machine)
JNIEnv *env(0);								/// Pointer to native method interface
jint JavaVersion = JNI_VERSION_1_6;			/// Minimum Java version required
#define CLEAR(x) memset(&x, 0, sizeof(x))	/// Clearer method


/// ************ PH_JAVAPLEX CLASS **********************************

std::pair<jobject, jobject> PH_Javaplex::computePH(Graph &G, jobject &stream)
{
	/// START CLOCK
	auto startTime = Utils::startChrono<std::chrono::high_resolution_clock::time_point>();
		
	/// Get maximum dimension, as Default(2) or BiggestClique (i.e. last clique in asc-sorted G.cliqueS)
	int cliquesNumber = (int) G.maximalCliqueSet.size();
	int maxDim = (cliquesNumber>0)? (int)G.maximalCliqueSet[cliquesNumber-1].size()+1 : 2+1;
	
	/// Set prime coefficient (2: Z2)
	int primeCoefficient = 2;
	
	/// Compute PH
	jobject persistence = getModularSimplicialAlgorithm(maxDim, primeCoefficient);
	//jobject persistence = Plex4_getModularSimplicialAlgorithm(plex4);

	/// Compute intervals
	jobject barcodes = computeIntervals(persistence, stream);	

	/// Compute intervals and generators (cycles)
	jobject annotatedBarcodes = computeAnnotatedIntervals(persistence, stream);
	//std::string barcode_string = getAnnotatedIntervals_AsString(annotatedBarcodes);
	//std::cout << "Raw barcodes: " << barcode_string << std::endl;

	/// Compute infinite intervals and generators (cycles)
	jobject infiniteBarcodes = getInfiniteIntervals(annotatedBarcodes);

	/// END CLOCK
	if(verbosity>0) std::cout << "Time(computePH): " << Utils::stopChrono(startTime) << "s" << std::endl;

	return std::make_pair(infiniteBarcodes, annotatedBarcodes);
}

jobject PH_Javaplex::j_buildComplex(Graph &G, int complexType, int dimension)
{
	/// Create simplex stream
	jobject stream = createExplicitSimplexStream();
	//jobject plex4 = Plex4();
	//jobject stream = Plex4_createExplicitSimplexStream(plex4);

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
			j_buildAllCliquesComplex(G, stream, autoAlgorithm, dimension);

			break;
		}
		case 1:					
		{
			/// Add Graph to stream
			/// Persistent homology: return <Elementary cycles> (AKA Holes AKA Hamiltonian) => minimal subgraph such that every node has degree 2
			///TODO: For G_008_complex some cycles are only expanded (e.g. 2 3 6 + 0 1; 1 5 7 + 0 2)
			j_buildGraphComplex(G, stream);

			break;
		}
		case 2:	
		{
			/// Add Maximal Cliques to stream
			/// Persistent homology: return <Elementary cycles>
			///TODO: For G_008_complex & G_012_cliques results are the same (9 & 15) as buildGraphComplex()
			j_buildMaxCliqueComplex(G, stream);
			
			break;
		}
		case 3:	
		{
			/// Add Minimum Cliques to stream
			/// Persistent homology: return <?>
			j_buildMinCyclesComplex(G, stream);
			
			break;
		}
	}

	/// Check consistency
	if(validateVerbose(stream))
	{
		if(verbosity>0) std::cout << "The complex is valid." << std::endl;
	}
	else
	{
		if(verbosity>0) std::cout << "The complex is not valid." << std::endl;
		ensureAllFaces(stream);
		//bool removed = removeElementIfPresent(stream, G.cliqueS[0]); //(dim+1 = G.V?)
	}

	/// Finalize stream
	finalizeStream(stream);

	/// Show filtration info
	/*if(verbosity > 0)
	{
		int streamsize = getSize(stream);
		std::cout << "\t*N. of vertices = " << G.V << std::endl;
		std::cout << "\t*N. of simplices = " << streamsize << std::endl;
		if(G.maximalCliqueSet.size() > 0)
			std::cout << "\t*Dimension = " << G.maximalCliqueSet[G.maximalCliqueSet.size()-1].size() << std::endl;
	}*/

	return stream;
}

void PH_Javaplex::j_buildGraphComplex(Graph &G, jobject &stream)
{
	/// Add vertex:		0-simplices
	for(int v = 0; v < G.V; ++v)
	{
		addVertex(stream, v);
		addElement(stream, std::vector<int>{v});
	}

	/// Add edges:		1-simplices
	for(int i = 0; i < G.edgesReference.size(); ++i)
	{
		int u = std::get<0>(G.edgesReference[i]);
		int v = std::get<1>(G.edgesReference[i]);
		//double w = std::get<2>(G.edgesReference[i]);
		//int edgeIndex = i;	/// "i" because edges are taken in order (nondecreasing)
		int filter = G.edgeRank(u, v);
		addElement(stream, std::vector<int>{u,v}, filter);
	}
}

void PH_Javaplex::j_buildMaxCliqueComplex(Graph &G, jobject &stream)
{
	/// Clique algorithm (k>=1) is chosen based on graph density
	if (G.density() >= 0.5)
		BronKerbosch(G, 1, 1);	///1:Tomita (better for dense graphs)
	else
		BronKerbosch(G, 2, 1);	///2:Epstein (better for sparse graphs)

	/// Add vertices
	for(int v = 0; v < G.V; ++v)
	{
		addVertex(stream, v);
	}

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
			addElement(stream, std::vector<int>{v});
		}
		else if(k == 2)
		{
			/// Add 1-simplices
			int u = clique.front();
			int v = clique.back();
			addElement(stream, std::vector<int>{u, v}, G.edgeRank(u, v));
		}
		else
		{
			for(int i = 0; i < k; ++i)
			{
				int u = clique[i];
				for(int j = i + 1; j < k; ++j)
				{
					int v = clique[j];
					/// Add 1-simplices
					addElement(stream, std::vector<int>{u, v}, G.edgeRank(u, v));		
				}
			}
		}
	}
}

void PH_Javaplex::j_buildMinCyclesComplex(Graph &G, jobject &stream)
{
	/// Compute MCB
	Amaldi(G);

	/// Add vertices
	for(int v = 0; v < G.V; ++v)
	{
		addVertex(stream, v);
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
		addElement(stream, cycleInt);
	}
}

void PH_Javaplex::j_buildAllCliquesComplex(Graph &G, jobject &stream, bool autoAlgorithm, int dimension)
{
	/// AllCliquesComplex(G);
	Homology::buildCliquesComplex(G, autoAlgorithm, dimension);

	/// Add filtration to simplicial complex
	//if(verbosity>0) std::cout << "Filtrations:" << std::endl;
	for (auto it = G.PH_Filtration.begin(); it != G.PH_Filtration.end(); ++it)
	{
		const int &rank = it->first;
		const auto &rankCliques = (it->second);
		//if(verbosity>0) std::cout << "\tFilter<" << rank << ">: ";
		bool firstPair = true;
		for (auto cliqueP : rankCliques)
		{
			std::vector<int> &clique = cliqueP.first;
			double weight = cliqueP.second;

			int k = (int)clique.size();

				if (k == 1)
				{
					/// Add 0-simplex
					addVertex(stream, clique[0]);
					/*if(verbosity>0)
					{
						if (!firstPair) std::cout << " ; ";
						std::cout << "[" << clique[0] << "]";
					}*/
				}
				else
				{
					/// Add (k-1)-simplex
					///OR ERROR HERE?
					addElement(stream, clique, rank);
					/*if(verbosity>0)
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
					}*/
				}
			firstPair = false;
		}
		//if(verbosity>0) std::cout << std::endl;
	}
}

std::vector<int> PH_Javaplex::getBettiNumbers(std::pair<jobject, jobject> &PH, bool onlyInfinite, bool getSequence)
{
	jobject &infiniteBarcodes = PH.first;
	jobject &annotatedBarcodes = PH.second;

	/// Betti numbers
	std::string bettiNumbers = (!onlyInfinite)? getBettiNumbers(annotatedBarcodes) : getBettiNumbers(infiniteBarcodes);
	std::cout << "Betti numbers:\t" << bettiNumbers << std::endl;

	/// Betti sequence
	std::vector<int> bettiSequence;
	if(getSequence)
	{
		bettiSequence = (!onlyInfinite)? getBettiSequence(annotatedBarcodes) : getBettiSequence(infiniteBarcodes);
		std::cout << "Betti sequence:\t";
		for(int i = 0; i < bettiSequence.size(); ++i)
			std::cout << bettiSequence[i] << " ";
		std::cout << std::endl;
	}
	
	return bettiSequence;
}

void PH_Javaplex::getIntervalsAndGenerators(Graph &G, std::pair<jobject, jobject> &PH, int dimension, int displayMode, bool onlyInfinite, bool saveCycles)
{
	jobject &infiniteBarcodes = PH.first;
	jobject &annotatedBarcodes = PH.second;

	std::vector<int> dimensionS = (dimension>=0)? std::vector<int>{dimension} : (onlyInfinite? getDimensions(infiniteBarcodes) 
																							 : getDimensions(annotatedBarcodes));

	/// Select mode for getting PH-Data
	switch (displayMode) {

		/// For 1-dim: Get <endpoints, generators, weight> && Add <cycles> (inf-holes) into graph

		/// Show <intervals,generators> multiline
		case 0:					
		{
			std::cout << "Intervals & Generators:" << std::endl;	
			for(std::size_t i = 0; i < dimensionS.size(); ++i)
			{
				std::multimap<std::string, std::string> IG_dim = (onlyInfinite)? getIntervalGeneratorPairsAtDimension(infiniteBarcodes, dimensionS[i]) 
																			: getIntervalGeneratorPairsAtDimension(annotatedBarcodes, dimensionS[i]);
				std::cout << "\tdim-" << dimensionS[i] << ": " << std::endl;
				for ( auto it = IG_dim.begin(); it != IG_dim.end(); ++it )
					std::cout << "\t" << it->first << " : " << it->second << std::endl;
			}

			break;
		}

		/// Get <intervals,generators> (in-one-line)
		case 1:
		{
			std::cout << "Intervals & Generators:" << std::endl;	
			for(std::size_t i = 0; i < dimensionS.size(); ++i)
			{
				std::string interval = (onlyInfinite)? getIntervalGeneratorPairsAtDimension_AsString(infiniteBarcodes, dimensionS[i]) 
										: getIntervalGeneratorPairsAtDimension_AsString(annotatedBarcodes, dimensionS[i]);
				std::cout << "\tdim-" << dimensionS[i] << ": " << interval << std::endl;
			}

			break;
		}
	}

	/// If specified, save cycles(infinite holes) into graph
	if(saveCycles)
	{
		barcodesIntoGraph(G, annotatedBarcodes);
		holesIntoGraph(G, infiniteBarcodes);
	}
}

void PH_Javaplex::holesIntoGraph(Graph &G, jobject infiniteBarcodes)
{
	/// Set dimension from where to get cycles
	int dim = 1;
	
	/// Get intervals from barcodes
	jobject intervals = getIntervalsAtDimension(infiniteBarcodes, dim);

	/// Get EndPoints from intervals
	bool skipInfinite = false;
	std::vector<std::pair<double, double>> endPoints = getEndPoints(intervals, dim, skipInfinite);			/// Unsorted
	//std::vector<std::pair<double, double>> endPoints = getEndPoints(infiniteBarcodes, dim, skipInfinite);	/// Sorted

	/// Get Generators from barcodes
	std::vector<std::vector<std::vector<int>>> G_dim = getGeneratorsAtDimension(infiniteBarcodes, dim);

	/// Graph PH storage as Multimap<Interval,Generators>
	std::multimap<std::pair<double,double>, std::vector<std::pair<int, int>>> PH_barcodes;

	/// Get <endpoints, generators, weight>
	if (verbosity > 1 && G_dim.size() > 0) std::cout << "Save holes@dim-1 into Graph: " << std::endl;
	for(std::size_t i = 0; i < G_dim.size(); ++i)
	{
		///Load generator into graph's MCB #Step1
		std::vector<ut> hole(G.E);
		double cycle_weight = 0;

		double tw = 0;
		if (verbosity > 1) std::cout << "\t -> Added cycle: ";

		///TODO: barcodes can also be saved from here (or from there). Cycles <=> endPoints[i][1] < inf
		/// interval "i"
		/*double inf = std::numeric_limits<double>::infinity(); ///999999
		if (verbosity > 1)
		{
			if(endPoints[i].first < inf)
				std::cout << "[";
			else
				std::cout << "(";
			std::cout << endPoints[i].first << ", " << endPoints[i].second;
			if(endPoints[i].second < inf)
				std::cout << "] : ";
			else
				std::cout << ") : ";
		}*/

		/// generator "i"
		for(std::size_t j = 0; j < G_dim[i].size(); ++j)
		{
			///TODO: more than 2 items
			int u = G_dim[i][j][0];
			int v = G_dim[i][j][1];
			if (verbosity > 1) std::cout << "(" << u << "," << v; 
			double w = G.edgeWeight(u, v);
			tw += w;
			//std::cout << ", w:" << w;
			if (verbosity > 1) std::cout << ") ";

			///Load generator into graph's MCB #Step2
			hole[G.edgeReference(u,v)] = 1;
		}

		if (verbosity > 1) std::cout << "\t w= " << tw << std::endl;

		///Load generator into graph's MCB #Step3
		if(hole.size()>2) G.addHoleToMCB(hole, cycle_weight);
	}
}

void PH_Javaplex::barcodesIntoGraph(Graph &G, jobject annotatedBarcodes)
{
	/// Get all dimensions having barcodes
	std::vector<int> dimensions = getDimensions(annotatedBarcodes);

	for(int d = 0; d < dimensions.size(); ++d)
	{
		/// Set dimension from where to get barcodes
		int dim = dimensions[d];
		DataStructures::barcodeCollection barcoll;

		/// Get intervals from barcodes
		jobject intervals = getIntervalsAtDimension(annotatedBarcodes, dim);

		/// Get EndPoints from intervals
		bool skipInfinite = false;
		std::vector<std::pair<double, double>> endPoints = getEndPoints(intervals, dim, skipInfinite);			/// Unsorted
		//std::vector<std::pair<double, double>> endPoints = getEndPoints(infiniteBarcodes, dim, skipInfinite);	/// Sorted

		/// Get Generators from barcodes
		std::vector<std::vector<std::vector<int>>> G_dim = getGeneratorsAtDimension(annotatedBarcodes, dim);

		/// Graph PH storage as Multimap<Interval,Generators>
		std::multimap<std::pair<double,double>, std::vector<std::pair<int, int>>> PH_barcodes;

		/// Get <endpoints, generators, weight>
		for(std::size_t i = 0; i < G_dim.size(); ++i)
		{
			/// Add to graph's persistent homology data
			barcoll.insert(std::make_pair(endPoints[i], G_dim[i]));
		}

		/// Add barcodes from current dimension
		G.PH_barcodes.insert(std::make_pair(dim, barcoll));
	}
}


/// ************ JavaPlex methods **********************************
jobject PH_Javaplex::createExplicitSimplexStream()
{
	jclass ExplicitSimplexStream = env->FindClass("edu/stanford/math/plex4/streams/impl/ExplicitSimplexStream");
	jmethodID ini_stream = env->GetMethodID(ExplicitSimplexStream, "<init>", "()V");		/// in: none; out: void
	
	///Create with NewObject()
	jobject stream = env->NewObject(ExplicitSimplexStream, ini_stream);	
	JVMinstance.verifyEXEC();

	/// Eventually, turn jclass into global reference
	//jobject globalRef_ExplicitSimplexStream = (jclass) env->NewGlobalRef(ExplicitSimplexStream);
	//env->DeleteLocalRef(ExplicitSimplexStream);

	return stream;
}

void PH_Javaplex::addVertex(jobject stream, const int &vertex, double filtrationValue)
{
	jclass ExplicitSimplexStream = env->FindClass("edu/stanford/math/plex4/streams/impl/ExplicitSimplexStream");
	JVMinstance.verifyEXEC();

	jmethodID addVertex;
	if(filtrationValue == 0)
		addVertex = env->GetMethodID(ExplicitSimplexStream, "addVertex", "(I)V");		/// in: integer; out: void
	else
		addVertex = env->GetMethodID(ExplicitSimplexStream, "addVertex", "(ID)V");		/// in: integer, double; out: void
	JVMinstance.verifyEXEC();

	if(filtrationValue == 0)
		env->CallVoidMethod(stream, addVertex, (jint)vertex);
	else
		env->CallVoidMethod(stream, addVertex, (jint)vertex, (jdouble)filtrationValue);
	JVMinstance.verifyEXEC();

	/*if(verbosity > 0)
		std::cout << "\tAdded vertex " << vertex << "(filter=" << filtrationValue << ")" << std::endl;*/
}

void PH_Javaplex::addElement(jobject stream, const std::vector<int> &basisElements, double filtrationValue)
{
	jclass ExplicitSimplexStream = env->FindClass("edu/stanford/math/plex4/streams/impl/ExplicitSimplexStream");
	JVMinstance.verifyEXEC();

	jmethodID addElement;
	if(filtrationValue == 0)
		addElement = env->GetMethodID(ExplicitSimplexStream, "addElement", "([I)V");	/// in: array of integers; out: void
	else
		addElement = env->GetMethodID(ExplicitSimplexStream, "addElement", "([ID)V");	/// in: array of integers, double; out: void
	JVMinstance.verifyEXEC();

	/// std::vector --> jintArray
	jintArray newarr = convert_intVector2jintArray(basisElements);

	/// double --> jdouble
	jdouble filter = (jdouble) filtrationValue;

	/// Add elems[] into stream
	/// NB: Duplicated elements are automatically ignored.
	if(filtrationValue == 0)
	{
		env->CallVoidMethod(stream, addElement, newarr);
	}
	else
	{
		env->CallVoidMethod(stream, addElement, newarr, filter);
	}
	JVMinstance.verifyEXEC();

	/// release the jintArray
	env->DeleteLocalRef(newarr);															
	JVMinstance.verifyEXEC();

	/*if(verbosity > 0)
	{
		std::cout << "\tAdded " << (nElems-1) << "-simplex [";
		for(auto &vertex : basisElements)
		{
			std::cout << vertex;
			if(vertex != basisElements.back())
				std::cout << ",";
		}
		std::cout << "]" << "(filter=" << filtrationValue << ")" << std::endl;
	}*/
}

void PH_Javaplex::finalizeStream(jobject stream)
{
	jclass ExplicitSimplexStream = env->FindClass("edu/stanford/math/plex4/streams/impl/ExplicitSimplexStream");
	JVMinstance.verifyEXEC();

	jmethodID finalizeStream = env->GetMethodID(ExplicitSimplexStream, "finalizeStream", "()V");	/// in: none; out: none
	JVMinstance.verifyEXEC();

	env->CallVoidMethod(stream, finalizeStream);
	JVMinstance.verifyEXEC();
}

int PH_Javaplex::getSize(jobject stream)
{
	jclass ExplicitSimplexStream = env->FindClass("edu/stanford/math/plex4/streams/impl/ExplicitSimplexStream");
	JVMinstance.verifyEXEC();

	jmethodID getSize = env->GetMethodID(ExplicitSimplexStream, "getSize", "()I");	/// in: none; out: int
	int i = env->CallIntMethod(stream, getSize);
	JVMinstance.verifyEXEC();

	return i;
}

jobject PH_Javaplex::getModularSimplicialAlgorithm(int maxDimension, int primeCoefficient)
{
	/// maxDimension: dimension in which homology will be computed (e.g. 3 means 0, 1, and 2)
	/// primeCoefficient: coefficient (any prime number) with which homology will be computed (e.g. 2 means Z=2Z)

	///ModularIntField.getInstance(prime)
	jclass ModularIntField = env->FindClass("edu/stanford/math/primitivelib/algebraic/impl/ModularIntField");
	JVMinstance.verifyEXEC();
	//jmethodID mintGetInstance = env->GetMethodID(ModularIntField, "<init>", "(I)V");
	//jobject mintField = env->NewObject(ModularIntField, mintGetInstance, (jint)primeCoefficient);
	jmethodID mintGetInstance = env->GetStaticMethodID(ModularIntField, "getInstance", "(I)Ledu/stanford/math/primitivelib/algebraic/impl/ModularIntField;");
	JVMinstance.verifyEXEC();
	jobject mintField = env->CallObjectMethod(ModularIntField, mintGetInstance, (jint)primeCoefficient);
	JVMinstance.verifyEXEC();

	///SimplexComparator.getInstance()
	jclass SimplexComparator = env->FindClass("edu/stanford/math/plex4/homology/chain_basis/SimplexComparator");
	JVMinstance.verifyEXEC();
	//jmethodID simplexCompGetInstance = env->GetMethodID(SimplexComparator, "<init>", "()V");
	//jobject simplexComp  = env->NewObject(SimplexComparator, simplexCompGetInstance);
	jmethodID simplexCompGetInstance = env->GetStaticMethodID(SimplexComparator, "getInstance", "()Ledu/stanford/math/plex4/homology/chain_basis/SimplexComparator;");
	JVMinstance.verifyEXEC();
	jobject simplexComp = env->CallObjectMethod(SimplexComparator, simplexCompGetInstance);
	JVMinstance.verifyEXEC();

	///IntAbsoluteHomology class
	jclass IntAbsoluteHomology = env->FindClass("edu/stanford/math/plex4/autogen/homology/IntAbsoluteHomology");
	JVMinstance.verifyEXEC();

	///IntAbsoluteHomology initializer
	jmethodID ini_persistence = env->GetMethodID(IntAbsoluteHomology, "<init>", "(Ledu/stanford/math/primitivelib/autogen/algebraic/IntAbstractField;Ljava/util/Comparator;II)V");
	JVMinstance.verifyEXEC();

	///(IntAbsoluteHomology) PersitenceHomology
	jobject persistence = env->NewObject(IntAbsoluteHomology, ini_persistence, mintField, simplexComp, (jint)0, (jint)maxDimension);
	JVMinstance.verifyEXEC();

	return persistence;
}

jobject PH_Javaplex::computeIntervals(jobject persistence, jobject stream)
{
	jclass AbstractPersistenceAlgorithm = env->FindClass("edu/stanford/math/plex4/homology/interfaces/AbstractPersistenceAlgorithm");
	JVMinstance.verifyEXEC();

	///public edu.stanford.math.plex4.homology.barcodes.BarcodeCollection<java.lang.Double> computeIntervals(edu.stanford.math.plex4.streams.interfaces.AbstractFilteredStream<T>);
	///descriptor: (Ledu/stanford/math/plex4/streams/interfaces/AbstractFilteredStream;)Ledu/stanford/math/plex4/homology/barcodes/BarcodeCollection;
	jmethodID computeIntervals = env->GetMethodID(AbstractPersistenceAlgorithm, "computeIntervals", "(Ledu/stanford/math/plex4/streams/interfaces/AbstractFilteredStream;)Ledu/stanford/math/plex4/homology/barcodes/BarcodeCollection;");
	JVMinstance.verifyEXEC();

	///circle_intervals = persistence.computeIntervals(stream)
	jobject bettiBarcodes = env->CallObjectMethod(persistence, computeIntervals, stream);
	JVMinstance.verifyEXEC();

	return bettiBarcodes;
}

jobject PH_Javaplex::getInfiniteIntervals(jobject barcodes)
{
	jclass AbstractPersistenceBasisAlgorithm = env->FindClass("edu/stanford/math/plex4/homology/barcodes/AnnotatedBarcodeCollection");
	JVMinstance.verifyEXEC();
	
	jmethodID computeAnnotatedIntervals = env->GetMethodID(AbstractPersistenceBasisAlgorithm, "getInfiniteIntervals", "()Ledu/stanford/math/plex4/homology/barcodes/AnnotatedBarcodeCollection;");
	JVMinstance.verifyEXEC();

	jobject intervalCycles = env->CallObjectMethod(barcodes, computeAnnotatedIntervals);
	JVMinstance.verifyEXEC();

	return intervalCycles;
}

std::string PH_Javaplex::getBettiNumbers(jobject barcodes)
{
	jclass BarcodeCollection = env->GetObjectClass(barcodes);
	JVMinstance.verifyEXEC();

	jmethodID getBettiNumbers = env->GetMethodID(BarcodeCollection, "getBettiNumbers", "()Ljava/lang/String;");
	JVMinstance.verifyEXEC();  

	jstring bettiNumbers = (jstring) env->CallObjectMethod(barcodes, getBettiNumbers);
	JVMinstance.verifyEXEC(); 
	
	return convert_jstring2string(bettiNumbers);
}

std::vector<int> PH_Javaplex::getBettiSequence(jobject barcodes)
{
	jclass BarcodeCollection = env->GetObjectClass(barcodes);
	JVMinstance.verifyEXEC();

	jmethodID getBettiSequence = env->GetMethodID(BarcodeCollection, "getBettiSequence", "()[I");
	JVMinstance.verifyEXEC();  

	jintArray bettiSequence = (jintArray) env->CallObjectMethod(barcodes, getBettiSequence);
	JVMinstance.verifyEXEC();
	
	return convert_jintArray2intVector(bettiSequence);
}

jobject PH_Javaplex::computeAnnotatedIntervals(jobject persistence, jobject stream)
{
	jclass AbstractPersistenceBasisAlgorithm = env->FindClass("edu/stanford/math/plex4/homology/interfaces/AbstractPersistenceBasisAlgorithm");
	JVMinstance.verifyEXEC();

	/// public edu.stanford.math.plex4.homology.barcodes.AnnotatedBarcodeCollection<java.lang.Double, B> computeAnnotatedIntervals(edu.stanford.math.plex4.streams.interfaces.AbstractFilteredStream<T>);
	/// descriptor: (Ledu/stanford/math/plex4/streams/interfaces/AbstractFilteredStream;)Ledu/stanford/math/plex4/homology/barcodes/AnnotatedBarcodeCollection;
	jmethodID computeAnnotatedIntervals = env->GetMethodID(AbstractPersistenceBasisAlgorithm, "computeAnnotatedIntervals", "(Ledu/stanford/math/plex4/streams/interfaces/AbstractFilteredStream;)Ledu/stanford/math/plex4/homology/barcodes/AnnotatedBarcodeCollection;");
	JVMinstance.verifyEXEC();

	///circle_intervals = persistence.computeAnnotatedIntervals(stream)
	jobject intervalCycles = env->CallObjectMethod(persistence, computeAnnotatedIntervals, stream);
	JVMinstance.verifyEXEC();

	return intervalCycles;
}

std::string PH_Javaplex::getAnnotatedIntervals_AsString(jobject annotatedBarcodes)
{
	jmethodID toString = env->GetMethodID(env->GetObjectClass(annotatedBarcodes), "toString", "()Ljava/lang/String;");
	JVMinstance.verifyEXEC();

	jstring result = (jstring) env->CallObjectMethod(annotatedBarcodes,toString);
	JVMinstance.verifyEXEC();

	return convert_jstring2string(result);
}

std::map<int, int> PH_Javaplex::getBettiNumbersMap(jobject annotatedBarcodes)
///BUG: not working
{
	jclass AnnotatedBarcodeCollection = env->GetObjectClass(annotatedBarcodes);
	JVMinstance.verifyEXEC();

	jmethodID getBettiNumbersMap = env->GetMethodID(AnnotatedBarcodeCollection, "getBettiNumbersMap", "(Ljava/lang/Comparable;)Ljava/util/Map;");
	JVMinstance.verifyEXEC();  

	jobject bettiNumbersMap = env->CallObjectMethod(annotatedBarcodes, getBettiNumbersMap, (jint)1);	/// param1: comparable
	JVMinstance.verifyEXEC(); 
	///Exception in thread "main" java.lang.NullPointerException
		///at edu.stanford.math.plex4.homology.barcodes.Interval.containsPoint(Unknown Source)
		///at edu.stanford.math.plex4.homology.barcodes.AnnotatedBarcodeCollection.getBettiNumbersMap(Unknown Source)

	jclass c_Map = env->FindClass("java/util/Map");
	jmethodID m_GetSize = env->GetMethodID(c_Map, "size", "()I");
	int jSize = env->CallIntMethod(bettiNumbersMap, m_GetSize);

	std::map<int, int> result;
	//env->ReleaseStringUTFChars(bettiNumbersMap, result_js);
	
	return result;
}

std::vector<int> PH_Javaplex::getDimensions(jobject annotatedBarcodes)
{
	jclass PersistenceInvariantDescriptor = env->FindClass("edu/stanford/math/plex4/homology/barcodes/PersistenceInvariantDescriptor");
	JVMinstance.verifyEXEC();

	jmethodID getDimensions = env->GetMethodID(PersistenceInvariantDescriptor, "getDimensions", "()Ljava/util/Set;");
	JVMinstance.verifyEXEC();

	jobject dimensions = env->CallObjectMethod(annotatedBarcodes, getDimensions);
	JVMinstance.verifyEXEC();

	jclass c_Set = env->GetObjectClass(dimensions);
	JVMinstance.verifyEXEC();

	jmethodID toArray = env->GetMethodID(c_Set, "toArray", "()[Ljava/lang/Object;");
	JVMinstance.verifyEXEC();

	jobjectArray arrayOfElements = (jobjectArray) env->CallObjectMethod(dimensions, toArray);
	JVMinstance.verifyEXEC();

	int arraySize = env->GetArrayLength(arrayOfElements);

	std::vector<int> result(arraySize);

	for (int i=0; i < arraySize; i++)
	{
		jobject obj = env->GetObjectArrayElement(arrayOfElements, i);
		jclass c_Integer = env->FindClass("java/lang/Integer");
		jmethodID intValue = env->GetMethodID(c_Integer, "intValue", "()I");
		result[i] = env->CallIntMethod(obj, intValue);;
	}

	return result;
}

std::string PH_Javaplex::getIntervalGeneratorPairsAtDimension_AsString(jobject annotatedBarcodes, int dim)
{
	jclass PersistenceInvariantDescriptor = env->FindClass("edu/stanford/math/plex4/homology/barcodes/PersistenceInvariantDescriptor");
	JVMinstance.verifyEXEC();

	jmethodID getIntervalGeneratorPairsAtDimension = env->GetMethodID(PersistenceInvariantDescriptor, "getIntervalGeneratorPairsAtDimension", "(I)Ljava/util/List;");
	JVMinstance.verifyEXEC();

	jobject intervalPairs= env->CallObjectMethod(annotatedBarcodes, getIntervalGeneratorPairsAtDimension, (jint)dim);
	JVMinstance.verifyEXEC();

	jmethodID ObjectObjectPair = env->GetMethodID(env->GetObjectClass(intervalPairs), "toString", "()Ljava/lang/String;");
	JVMinstance.verifyEXEC();

	jstring result = (jstring) env->CallObjectMethod(intervalPairs,ObjectObjectPair);
	JVMinstance.verifyEXEC();
	
	return convert_jstring2string(result);
}

std::multimap<std::string, std::string> PH_Javaplex::getIntervalGeneratorPairsAtDimension(jobject annotatedBarcodes, int dim)
{
	jclass PersistenceInvariantDescriptor = env->FindClass("edu/stanford/math/plex4/homology/barcodes/PersistenceInvariantDescriptor");
	JVMinstance.verifyEXEC();

	jmethodID getIntervalGeneratorPairsAtDimension = env->GetMethodID(PersistenceInvariantDescriptor, "getIntervalGeneratorPairsAtDimension", "(I)Ljava/util/List;");
	JVMinstance.verifyEXEC();

	///Return: List<ObjectObjectPair<I, G>>
	///https://github.com/ogdet/primitive-lib/blob/master/src/edu/stanford/math/primitivelib/autogen/pair/ObjectObjectPair.java
	jobject intervalPairs = env->CallObjectMethod(annotatedBarcodes, getIntervalGeneratorPairsAtDimension, (jint)dim);
	JVMinstance.verifyEXEC();

	jclass List = env->GetObjectClass(intervalPairs);
	JVMinstance.verifyEXEC();

	///List.size
	jmethodID size = env->GetMethodID(List, "size", "()I");
	JVMinstance.verifyEXEC();
	int n = env->CallIntMethod(intervalPairs, size);
	JVMinstance.verifyEXEC();

	///OOPair getFirst & getSecond
	jclass ObjectObjectPair = env->FindClass("edu/stanford/math/primitivelib/autogen/pair/ObjectObjectPair");
	JVMinstance.verifyEXEC();
	jmethodID getFirst = env->GetMethodID(ObjectObjectPair, "getFirst", "()Ljava/lang/Object;");
	JVMinstance.verifyEXEC();
	jmethodID getSecond = env->GetMethodID(ObjectObjectPair, "getSecond", "()Ljava/lang/Object;");
	JVMinstance.verifyEXEC();
	///List.get
	jmethodID get = env->GetMethodID(List, "get", "(I)Ljava/lang/Object;" );
	JVMinstance.verifyEXEC();

	///Iterate List items
	std::multimap<std::string, std::string> result;
	for( int i = 0; i < n; ++i )
	{
		///Return: ObjectObjectPair<I, G>
		jobject item = env->CallObjectMethod(intervalPairs, get, i );
		JVMinstance.verifyEXEC();
		
		///Return: I_nterval = HashMap<Integer, List<I> e.g. [0,inf)
		jstring item_I = (jstring) env->CallObjectMethod(item, getFirst);
		JVMinstance.verifyEXEC();
		jclass I_class = env->GetObjectClass(item_I);
		JVMinstance.verifyEXEC();
		jmethodID ItoString = env->GetMethodID(I_class, "toString", "()Ljava/lang/String;");
		JVMinstance.verifyEXEC();
		jstring I_jstring = (jstring) env->CallObjectMethod(item_I, ItoString);
		JVMinstance.verifyEXEC();

		///Return: G_enerators = HashMap<Integer, List<G> e.g. [4,5] + [0,4] + [0,5]
		jobject item_G = (jobject) env->CallObjectMethod(item, getSecond);
		JVMinstance.verifyEXEC();
		jclass G_class = env->GetObjectClass(item_G);
		JVMinstance.verifyEXEC();
		jmethodID GtoString = env->GetMethodID(G_class, "toString", "()Ljava/lang/String;");
		JVMinstance.verifyEXEC();
		jstring G_jstring = (jstring) env->CallObjectMethod(item_G, GtoString);
		JVMinstance.verifyEXEC();
		
		result.insert(std::make_pair(convert_jstring2string(I_jstring), convert_jstring2string(G_jstring)));
	}
		
	return result;
}

jobject PH_Javaplex::getIntervalsAtDimension(jobject annotatedBarcodes, int dim)
{
	jclass PersistenceInvariantDescriptor = env->FindClass("edu/stanford/math/plex4/homology/barcodes/PersistenceInvariantDescriptor");
	JVMinstance.verifyEXEC();

	///public java.util.List<I> getIntervalsAtDimension(int dimension)
	jmethodID getIntervalsAtDimension = env->GetMethodID(PersistenceInvariantDescriptor, "getIntervalsAtDimension", "(I)Ljava/util/List;");
	JVMinstance.verifyEXEC();  

	jobject intervals = env->CallObjectMethod(annotatedBarcodes, getIntervalsAtDimension, (jint)dim);
	JVMinstance.verifyEXEC();

	return intervals;
}

std::vector<std::vector<std::vector<int>>> PH_Javaplex::getGeneratorsAtDimension(jobject annotatedBarcodes, int dim)
{
	///List<G> getGeneratorsAtDimension(int dimension)
	jclass PersistenceInvariantDescriptor = env->FindClass("edu/stanford/math/plex4/homology/barcodes/PersistenceInvariantDescriptor");
	JVMinstance.verifyEXEC();

	jmethodID getGeneratorsAtDimension = env->GetMethodID(PersistenceInvariantDescriptor, "getGeneratorsAtDimension", "(I)Ljava/util/List;");
	JVMinstance.verifyEXEC();

	jobject generators = env->CallObjectMethod(annotatedBarcodes, getGeneratorsAtDimension, (jint)dim);
	JVMinstance.verifyEXEC();

	jclass List = env->GetObjectClass(generators);
	JVMinstance.verifyEXEC();

	///List.size
	jmethodID size = env->GetMethodID(List, "size", "()I");
	JVMinstance.verifyEXEC();
	int n = env->CallIntMethod(generators, size);
	JVMinstance.verifyEXEC();

	///List.get
	jmethodID get = env->GetMethodID(List, "get", "(I)Ljava/lang/Object;" );
	JVMinstance.verifyEXEC();

	///Iterate List items
	std::vector<std::vector<std::vector<int>>> result;
	for( int i = 0; i < n; ++i )
	{
		///Return: dim-1 [4,5] + [0,4] + [0,5]
		///        dim-2 [4,5,8] + [2,7,10] + [5,8,13]
		jobject gen = env->CallObjectMethod(generators, get, i);
		JVMinstance.verifyEXEC();

		jclass c_gen = env->GetObjectClass(gen);
		JVMinstance.verifyEXEC();

		jmethodID toString = env->GetMethodID(List, "toString", "()Ljava/lang/String;");
		JVMinstance.verifyEXEC();

		jstring gen_jstring = (jstring) env->CallObjectMethod(gen, toString, i);
		JVMinstance.verifyEXEC();

		result.push_back(convert_generator2vector(convert_jstring2string(gen_jstring)));
	}

	return result;
}

std::vector<std::pair<double, double>> PH_Javaplex::getEndPoints(jobject annotatedBarcodes, int dim, bool skipInfIntervals)
{
	/// Get BarcodeUtility class
	jclass BarcodeUtility = env->FindClass("edu/stanford/math/plex4/homology/barcodes/BarcodeUtility");
	JVMinstance.verifyEXEC();

	/// Create new instance of BarcodeUtility
	jmethodID ini_barcodeUtility = env->GetMethodID(BarcodeUtility, "<init>", "()V");
	JVMinstance.verifyEXEC();
	jobject barcodeUtility = env->NewObject(BarcodeUtility, ini_barcodeUtility);
	JVMinstance.verifyEXEC();

	/// Get class of input param
	jclass myClass = env->GetObjectClass(annotatedBarcodes);
	JVMinstance.verifyEXEC();

	/// Get possible matching classes
	jclass annotatedBarcodes_class = env->GetObjectClass(annotatedBarcodes);
	JVMinstance.verifyEXEC();
	jclass intervalsList_class = env->FindClass("Ljava/util/List;");
	JVMinstance.verifyEXEC();

	/// Get method according to object class
	/// method: double[][] getEndpoints(AnnotatedBarcodeCollection<Double, G>, int, boolean)
	jmethodID getEndpoints;
	if(myClass == annotatedBarcodes_class)
		getEndpoints = env->GetStaticMethodID(BarcodeUtility, "getEndpoints", "(Ledu/stanford/math/plex4/homology/barcodes/AnnotatedBarcodeCollection;IZ)[[D");
	else
		getEndpoints = env->GetStaticMethodID(BarcodeUtility, "getEndpoints", "(Ljava/util/List;IZ)[[D");
	JVMinstance.verifyEXEC();
	
	jobject endPoints= env->CallObjectMethod(barcodeUtility, getEndpoints, annotatedBarcodes, (jint)dim, (jboolean)skipInfIntervals);
	JVMinstance.verifyEXEC();

	///Convert jdouble[][] to vector<vector<double>>
	jobjectArray arrayOfElements = (jobjectArray)endPoints;
	JVMinstance.verifyEXEC();
	int size = env->GetArrayLength(arrayOfElements);

	///Iterate List items
	//std::vector<std::vector<double>> result(size);
	std::vector<std::pair<double, double>> result(size);
	for( int i = 0; i < size; ++i )
	{
		jdoubleArray endPoint = (jdoubleArray)env->GetObjectArrayElement(arrayOfElements, i);
		std::vector<double> temp = convert_jdoubleArray2doubleVector(endPoint);
		std::pair<double, double> interval = std::make_pair(temp[0], temp[1]);
		//result[i] = temp;
		result[i] = interval;
	}
	
	return result;
}

void PH_Javaplex::ensureAllFaces(jobject stream)
{
	jclass ExplicitSimplexStream = env->FindClass("edu/stanford/math/plex4/streams/impl/ExplicitSimplexStream");
	JVMinstance.verifyEXEC();

	jmethodID ensureAllFaces = env->GetMethodID(ExplicitSimplexStream, "ensureAllFaces", "()V");
	JVMinstance.verifyEXEC();

	env->CallVoidMethod(stream, ensureAllFaces);
	JVMinstance.verifyEXEC();
}

bool PH_Javaplex::removeElementIfPresent(jobject stream, const std::vector<int> &basisElements)
{
	jclass ExplicitSimplexStream = env->FindClass("edu/stanford/math/plex4/streams/impl/ExplicitSimplexStream");
	JVMinstance.verifyEXEC();

	jmethodID removeElementIfPresent = env->GetMethodID(ExplicitSimplexStream, "removeElementIfPresent", "([I)Z");
	JVMinstance.verifyEXEC();

	/// std::vector --> jintArray
	jintArray newarr = convert_intVector2jintArray(basisElements);

	jboolean result = env->CallBooleanMethod(stream, removeElementIfPresent, newarr);
	JVMinstance.verifyEXEC();

	/// release the jintArray
	env->DeleteLocalRef(newarr);															
	JVMinstance.verifyEXEC();

	return result;
}

bool PH_Javaplex::validateVerbose(jobject stream)
{
	jclass ExplicitSimplexStream = env->FindClass("edu/stanford/math/plex4/streams/impl/ExplicitSimplexStream");
	JVMinstance.verifyEXEC();

	///Validates the stream to make sure that it contains a valid filtered simplicial complex. It checks:
	/// 1. For each element in the complex, all of the faces of the simplex also belong to the complex.
	/// 2. The faces of each simplex have filtration values that are less than or equal to those of its cofaces.
	jmethodID validateVerbose = env->GetMethodID(ExplicitSimplexStream, "validateVerbose", "()Z");
	JVMinstance.verifyEXEC();

	jboolean result = env->CallBooleanMethod(stream, validateVerbose);
	JVMinstance.verifyEXEC(); 

	return result;
}

std::map<std::vector<std::vector<int>>, double> PH_Javaplex::createBoundaryMatrixAsDoubleSum(jobject stream)
{
	/// Get StreamUtility class
	jclass StreamUtility = env->FindClass("edu/stanford/math/plex4/streams/utility/StreamUtility");
	JVMinstance.verifyEXEC();

	///public static <T> edu.stanford.math.primitivelib.autogen.formal_sum.DoubleSparseFormalSum<edu.stanford.math.primitivelib.autogen.pair.ObjectObjectPair<T, T>> createBoundaryMatrixAsDoubleSum(edu.stanford.math.plex4.streams.interfaces.AbstractFilteredStream<T>);
	///descriptor: (Ledu/stanford/math/plex4/streams/interfaces/AbstractFilteredStream;)Ledu/stanford/math/primitivelib/autogen/formal_sum/DoubleSparseFormalSum;	
	jmethodID createBoundaryMatrixAsDoubleSum = env->GetStaticMethodID(StreamUtility, "createBoundaryMatrixAsDoubleSum", "(Ledu/stanford/math/plex4/streams/interfaces/AbstractFilteredStream;)Ledu/stanford/math/primitivelib/autogen/formal_sum/DoubleSparseFormalSum;");
	JVMinstance.verifyEXEC();

	/// Get Boundary Matrix as "DoubleSparseFormalSum<ObjectObjectPair<T, T>>"
	jobject boundaryMatrix = env->CallObjectMethod(stream, createBoundaryMatrixAsDoubleSum, stream);
	JVMinstance.verifyEXEC();

	/// Get class ("edu/stanford/math/primitivelib/autogen/formal_sum/DoubleSparseFormalSum")
	///https://github.com/ogdet/primitive-lib/blob/master/src/edu/stanford/math/primitivelib/autogen/formal_sum/DoubleSparseFormalSum.java
	jclass BM_class = env->GetObjectClass(boundaryMatrix);
	JVMinstance.verifyEXEC();

	/// Convert to string (List <pair<Double coefficient,TObjectDoubleHashMap map>)
	jmethodID toString = env->GetMethodID(BM_class, "toString", "()Ljava/lang/String;");
	JVMinstance.verifyEXEC();
	jstring boundaryMatrix_jstring = (jstring) env->CallObjectMethod(boundaryMatrix, toString);
	JVMinstance.verifyEXEC();
	std::string boundaryMatrixS = convert_jstring2string(boundaryMatrix_jstring);

	/// Get size
	jmethodID size = env->GetMethodID(BM_class, "size", "()I");
	JVMinstance.verifyEXEC();
	int n = env->CallIntMethod(boundaryMatrix, size);
	JVMinstance.verifyEXEC();

	/// Get iterator (DoubleSparseFormalSum.map.iterator ;where map is TObjectDoubleHashMap<M>)
	jmethodID BM_iterator_method = env->GetMethodID(BM_class, "iterator", "()Lgnu/trove/TObjectDoubleIterator;");
	JVMinstance.verifyEXEC();
	jobject BM_iterator = env->CallObjectMethod(boundaryMatrix, BM_iterator_method);
	JVMinstance.verifyEXEC();
	jclass BM_iterator_class = env->GetObjectClass(BM_iterator);
	JVMinstance.verifyEXEC();

	/// Get iterator methods
	jmethodID hasNext = env->GetMethodID(BM_iterator_class, "hasNext", "()Z");
	JVMinstance.verifyEXEC();

	jmethodID advance = env->GetMethodID(BM_iterator_class, "advance", "()V");
	JVMinstance.verifyEXEC();
	
	jmethodID value = env->GetMethodID(BM_iterator_class, "value", "()D");
	JVMinstance.verifyEXEC();	

	jmethodID key = env->GetMethodID(BM_iterator_class, "key", "()Ljava/lang/Object;");
	JVMinstance.verifyEXEC();
		
	///Iterate List items
	std::map<std::vector<std::vector<int>>, double> result;
	while (env->CallBooleanMethod(BM_iterator, hasNext))
	{
		///Return:		 1 ([15] [7,15]
		///				-1 [1] [1,8]

		env->CallVoidMethod(BM_iterator, advance);
		JVMinstance.verifyEXEC();
		
		///Coefficient: Double
		jdouble item_value = env->CallDoubleMethod(BM_iterator, value);
		JVMinstance.verifyEXEC();

		///Object: ObjectObjectPair
		jobject item_key = env->CallObjectMethod(BM_iterator, key);
		JVMinstance.verifyEXEC();
			
		/// Get items separately. They do not have "size", "hasNext"
		/*///item_key_A
		jmethodID getFirst = env->GetMethodID(env->GetObjectClass(item_key), "getFirst", "()Ljava/lang/Object;");
		JVMinstance.verifyEXEC();
		jobject item_key_A = env->CallObjectMethod(item_key, getFirst);
		JVMinstance.verifyEXEC();
		///item_key_B
		jmethodID getSecond = env->GetMethodID(env->GetObjectClass(item_key), "getSecond", "()Ljava/lang/Object;");
		JVMinstance.verifyEXEC();
		jobject item_key_B = env->CallObjectMethod(item_key, getSecond);
		JVMinstance.verifyEXEC();*/

		///Get <Object,Coefficient> into result
		jmethodID itemtoString = env->GetMethodID(env->GetObjectClass(item_key), "toString", "()Ljava/lang/String;");
		JVMinstance.verifyEXEC();
		jstring item_jstring = (jstring) env->CallObjectMethod(item_key, itemtoString);
		JVMinstance.verifyEXEC();
		
		std::vector<std::vector<int>> _key = convert_generator2vector(convert_jstring2string(item_jstring));
		double _value = item_value;
		result.insert(std::make_pair(_key, _value));
	}	

	return result;
}

/// ************ PLEX4 CLASS **********************************

jobject PH_Javaplex::Plex4()
{
	/// Get Plex4 class
	jclass apiPlex4 = env->FindClass("edu/stanford/math/plex4/api/Plex4");
	JVMinstance.verifyEXEC();
	
	/// Create new instance of Plex4
	jmethodID ini_plex4 = env->GetMethodID(apiPlex4, "<init>", "()V");
	jobject plex4 = env->NewObject(apiPlex4, ini_plex4);
	JVMinstance.verifyEXEC();
	
	///Retain a reference to object for use in future (i.e. no need to pass it as a JNI method parameter). Remember: DeleteGlobalRef()
	env->NewGlobalRef(plex4);

	return plex4;
}

jobject PH_Javaplex::Plex4_createExplicitSimplexStream(jobject plex4)
{
	/// Get Plex4 class
	jclass apiPlex4 = env->FindClass("edu/stanford/math/plex4/api/Plex4");
	JVMinstance.verifyEXEC();

	/// Create simplex
	jmethodID createExplicitSimplexStream = env->GetStaticMethodID(apiPlex4, "createExplicitSimplexStream", "()Ledu/stanford/math/plex4/streams/impl/ExplicitSimplexStream;");
	JVMinstance.verifyEXEC();

	jobject stream = env->CallObjectMethod(plex4, createExplicitSimplexStream);
	JVMinstance.verifyEXEC();
	
	///Retain a reference to object for use in future (i.e. no need to pass it as a JNI method parameter). Remember: DeleteGlobalRef()
	env->NewGlobalRef(stream);

	return stream;
}

jobject PH_Javaplex::Plex4_getModularSimplicialAlgorithm(jobject plex4)
{
	/// Get Plex4 class
	jclass apiPlex4 = env->FindClass("edu/stanford/math/plex4/api/Plex4");
	JVMinstance.verifyEXEC();

	plex4 = env->NewGlobalRef(plex4);

	/// Get algorithm
	jmethodID getModularSimplicialAlgorithm = env->GetStaticMethodID(apiPlex4, "getModularSimplicialAlgorithm", "(II)Ledu/stanford/math/plex4/homology/interfaces/AbstractPersistenceAlgorithm;");
	JVMinstance.verifyEXEC();
	jobject persistence = env->CallObjectMethod(plex4, getModularSimplicialAlgorithm, (jint)3, (jint)2);
	JVMinstance.verifyEXEC();
	
	///Retain a reference to object for use in future (i.e. no need to pass it as a JNI method parameter). Remember: DeleteGlobalRef()
	env->NewGlobalRef(persistence);

	return persistence;
}

/// ************ Conversion methods **********************************

std::string PH_Javaplex::convert_jstring2string(jstring j_string)
{
	const char *c_chararray = env->GetStringUTFChars(j_string, NULL);
	std::string c_string(c_chararray );
	env->ReleaseStringUTFChars(j_string, c_chararray);

	return c_string;
}

jintArray PH_Javaplex::convert_intVector2jintArray(const std::vector<int> &basisElements)
{
	std::vector<jint> tmp(basisElements.begin(), basisElements.end());
	jint *arr = &tmp[0];
	const jsize nElems = (jsize) basisElements.size();
	jintArray newarr = env->NewIntArray(nElems);
	env->SetIntArrayRegion(newarr, 0, nElems, arr);
	JVMinstance.verifyEXEC();

	///env->DeleteLocalRef(newarr);	

	return newarr;
}

std::vector<int> PH_Javaplex::convert_jintArray2intVector(const jintArray &basisElements)
{
	const jsize length = env->GetArrayLength(basisElements);
	jint *bseq = env->GetIntArrayElements(basisElements, NULL);

	std::vector<int> result(length);

	for (int i = 0; i < length; i++)
	{
		result[i] = bseq[i];
	}

	env->ReleaseIntArrayElements(basisElements, bseq, NULL);

	return result;
}

jdoubleArray PH_Javaplex::convert_doubleVector2jdoubleArray(const std::vector<double> &basisElements)
{
	std::vector<jdouble> tmp(basisElements.begin(), basisElements.end());
	jdouble *arr = &tmp[0];
	const jsize nElems = (jsize) basisElements.size();
	jdoubleArray newarr = env->NewDoubleArray(nElems);
	env->SetDoubleArrayRegion(newarr, 0, nElems, arr);
	JVMinstance.verifyEXEC();

	///env->DeleteLocalRef(newarr);	

	return newarr;
}

std::vector<double> PH_Javaplex::convert_jdoubleArray2doubleVector(const jdoubleArray &basisElements)
{
	const jsize length = env->GetArrayLength(basisElements);
	jdouble *bseq = env->GetDoubleArrayElements(basisElements, NULL);

	std::vector<double> result(length);

	for (int i = 0; i < length; i++)
	{
		result[i] = bseq[i];
	}

	env->ReleaseDoubleArrayElements(basisElements, bseq, NULL);

	return result;
}

std::vector<std::vector<int>> PH_Javaplex::convert_generator2vector(std::string generatorS)
{
	/// INPUT: [4,5] + [0,4] + [0,5]
	/// OUTPUT: { (4,5), (0,4), (0,5) }

	/// INPUT: [4,5,7] + [0,4,9] + [1, 0,5]
	/// OUTPUT: { (4,5, 7), (0,4, 9), (1, 0,5) }

	std::vector<std::vector<int>> result;
	std::vector<int> vec;
	std::string number;

	/// Read each character in the string
	for(char& c : generatorS)
	{
		///element starts
		if(c == '[')
		{
			vec.clear();
		}
		else
		{
			if(c == ',')
			{
				vec.push_back(atoi(number.c_str()));
				number.clear();
			}
			else
			{
				if(c == ']')
				{
					vec.push_back(atoi(number.c_str()));
					number.clear();
					result.push_back(vec);
				}
				else
				{
					if(isdigit(c))
					{
						number.append(std::to_string((c-'0')));
					}
				}
			}
		}
	}

	return result;
}

/// ************ JVM CLASS **********************************
jsize JVM::GetnVMs()
{
	jvmHandle = GetModuleHandleA ("jvm.dll");
	if(jvmHandle == NULL)
	{
		Utils::writeLog("INFO. GetnVMs -> jvmHandle is null");
		return -1;
	}
	else
	{
		Utils::writeLog("INFO. GetnVMs -> jvmHandle is ok");
		FARPROC func_JNI_GetCreatedJavaVMs = GetProcAddress(jvmHandle , "JNI_GetCreatedJavaVMs");
		typedef jint (*hJNI_GetCreatedJavaVMs )( JavaVM** vmBuf , jsize bufLen , jsize* nVMs );		
		hJNI_GetCreatedJavaVMs oJNI_GetCreatedJavaVMs = (hJNI_GetCreatedJavaVMs) func_JNI_GetCreatedJavaVMs;
	
		jsize nVMs;
		oJNI_GetCreatedJavaVMs(NULL, 0, &nVMs);
		JavaVM** buffer = new JavaVM*[nVMs];
		oJNI_GetCreatedJavaVMs(buffer, nVMs, &nVMs);

		Utils::writeLog("INFO. GetnVMs -> " + (int)nVMs);
		return nVMs;
	}
}

JNIEnv* JVM::CreateJVM()
{ 
	/// PREPARE JVM initial options
	CLEAR(jvm);
	CLEAR(env);
	JavaVMOption options;
	CLEAR(options);
	std::string& option1 = std::string("-Djava.class.path=");				/// location of jar/class file(s)
	option1 += ".";															/// application directory (same as .\\)
	option1 += ";.\\javaplex.jar";											/// javaplex jar file
	//option1 += ";D:\\Graph\\graph-suite\\graph-app\\bin\\Release\\javaplex.jar";
	//const std::string& option2 = std::string("-verbose:jni");				/// print JNI-related messages
	//const std::string& option3 = std::string("-Djava.compiler=NONE");		/// disable JIT
	//const std::string& option4 = std::string("-Djava.library.path=.");		/// native library path	
	std::string optionString = option1;// + ";" + option2 + ";" + option3 + ";" + option4;
	options.optionString = const_cast<char*>( optionString.c_str() );
	options.extraInfo = 0;

	/// PREPARE JVM initial arguments
	JavaVMInitArgs vm_args;
	CLEAR(vm_args);
	memset( &vm_args, 0, sizeof( vm_args ) );		/// required to prevent memory issues
	vm_args.version = JavaVersion;					/// minimum java version
	vm_args.nOptions = 1;							/// number of options
	vm_args.options = &options;						/// set options
	vm_args.ignoreUnrecognized = JNI_TRUE;			/// JNI_TRUE: Skip wrong options; JNI_FALSE: Exit if wrong options

	Utils::writeLog("INFO. JVM arguments: " + static_cast<std::string>(vm_args.options->optionString));

	/// GET JVM.DLL PATH from H_key_LocalMAchine
	std::string javaVersionA = Utils::GetStringRegKey("SOFTWARE\\JavaSoft\\Java Runtime Environment\\", "CurrentVersion");
	std::string javaVersionB = Utils::GetStringRegKey("SOFTWARE\\JavaSoft\\JRE\\", "CurrentVersion");
	std::string javaPathA = Utils::GetStringRegKey("SOFTWARE\\JavaSoft\\Java Runtime Environment\\" + javaVersionA + "\\", "RuntimeLib");
	std::string javaPathB = Utils::GetStringRegKey("SOFTWARE\\JavaSoft\\JRE\\" + javaVersionB + "\\", "RuntimeLib");
	std::string javaPath = ( javaPathA != "")? javaPathA : javaPathB;
	std::string javaVersion = ( javaVersionA != "")? javaVersionA : javaVersionB;
	std::wstring tmp = Utils::ConvertStringToWString(javaPath);
	auto jvmFullPath = tmp.c_str();
	if(javaVersion == "" || javaVersion == "0")
		Utils::writeLog("ERROR. No version of Java was detected!");
	else
		Utils::writeLog("OK. Java detected: " + javaVersion + " at " + javaPath);
	
	///CHECK for running JVM instances
	jint result;
	int nVMs = GetnVMs();
	if(nVMs < 0)
	{
		/// LOAD JVM and INITIALIZE JNI
		//hinstLib = LoadLibrary(TEXT("C:\\Program Files\\Java\\jdk1.8.0_162\\jre\\bin\\server\\jvm.dll"));
		//hinstLib = LoadLibrary(TEXT("jvm.dll"));
		if (hinstLib == 0 || hinstLib == NULL)
			hinstLib = LoadLibrary(jvmFullPath);
		///TODO: Last bug here, from AWS_WIN2016
		if (hinstLib == 0 || hinstLib == NULL)		
			Utils::writeLog("ERROR. Either the library cannot be loaded or it's NULL.");
		else
			Utils::writeLog("OK. JVM library successfully loaded.");
		//Sleep(1000);	/// Give 1s to the library to get loaded	
	}
	if(nVMs <= 0)
	{
		///CreateJVM
		typedef jint(JNICALL *pCreateJavaVM)(JavaVM **, void**, void *);
		pCreateJavaVM CreateJavaVM = (pCreateJavaVM)GetProcAddress(hinstLib, "JNI_CreateJavaVM");
		result = CreateJavaVM(&jvm, (void**)&env, &vm_args);			/// Dynamic way
		//jint result = JNI_CreateJavaVM(&jvm, (void**)&env, &vm_args);	/// Static way
		/// SEGV (exception 0xC0000005) is generated intentionally on JVM startup to verify certain CPU/OS features.
		/// Disable: Exception Settings > Win32 Exceptions > C00000005 Access violation
	}
	else
	{
		result = 666;
	}

	/// CHECK JVM status: JNI_OK(0) JNI_ERR(-1) JNI_EDETACHED(-2) JNI_EVERSION(-3) JNI_ENOMEM(-4) JNI_EEXIST(-5) JNI_EINVAL(-6)
	if ( result >= 0)
	{
		/// DISPLAY version
		Utils::writeLog("OK. CreateJavaVM succeeded. Version: " + ((env->GetVersion()>>16)&0x0f) + '.' + (env->GetVersion()&0x0f));

		/// EXECUTE code
		//calculateHomology();
		
		/// UNLOAD JVM
		//FreeLibrary(hinstLib);
		//jvm->DestroyJavaVM();
		
		/// RETURN environment
		return env;
	}
	else
	{
		/// throw std::exception("Error.\n");
		if(result == JNI_EVERSION)
			Utils::writeLog("ERROR. CreateJavaVM -> JVM is outdated and doesn't meet requirements");
		else if(result == JNI_ENOMEM)
			Utils::writeLog("ERROR. CreateJavaVM -> not enough memory for JVM");
		else if(result == JNI_EINVAL)
			Utils::writeLog("ERROR. CreateJavaVM -> invalid argument for launching JVM");
		else if(result == JNI_EEXIST)
			Utils::writeLog("ERROR. CreateJavaVM -> the process can only launch one JVM instance");
		else if(result == 666)
			Utils::writeLog("INFO. CreateJavaVM -> creation skipped. A JVM instance is already running");
		else
			Utils::writeLog("ERROR. CreateJavaVM -> could not create the JVM instance (error code " + result + ')');
		return NULL;
	}
}

JNIEnv* JVM::AttachJVM()
{
	/// ATTACH to running JVM
	if (hinstLib == 0 || hinstLib == NULL) 
	{
		/// GET JVM.DLL PATH
		std::string javaVersion = Utils::GetStringRegKey("SOFTWARE\\JavaSoft\\Java Runtime Environment\\", "CurrentVersion");
		std::string javaPath = Utils::GetStringRegKey("SOFTWARE\\JavaSoft\\Java Runtime Environment\\" + javaVersion + "\\", "RuntimeLib");
		std::wstring tmp = Utils::ConvertStringToWString(javaPath);
		auto jvmFullPath = tmp.c_str();

		/// LOAD JVM
		//hinstLib = LoadLibrary(TEXT("C:\\Program Files\\Java\\jdk1.8.0_162\\jre\\bin\\server\\jvm.dll"));
		//hinstLib = LoadLibrary(TEXT("jvm.dll"));
		hinstLib = LoadLibrary(jvmFullPath);
	}

	jint _attachNew = jvm->AttachCurrentThread ( ( void ** ) &env , NULL );

	if (env != nullptr)
	{
		Utils::writeLog("OK. AttachJVM succeed.");
		/// EXECUTE code
		//calculateHomology();

		/// DETACH Java VM
		//jvm->DetachCurrentThread( );

		/// RETURN environment
		return env;
	}
	else
	{
		/// throw std::exception("Error: unable to attach to Java Virtual Machine.\n");
		Utils::writeLog("ERROR. unable to attach to Java Virtual Machine.");
		return NULL;
	}
}

void JVM::StartJVM()
{
	jvmHandle = GetModuleHandleA ("jvm.dll");
	if(jvmHandle == NULL)
	{
		Utils::writeLog("INFO. StartJVM -> CreateJVM starting...");
		CreateJVM();
	}
	else
	{
		/// CHECK for running JVM
		jsize nVMs = GetnVMs();

		/// RUNNING instance found
		if(nVMs > 0)
		{
			jint _env = jvm->GetEnv(( void ** ) &env, JavaVersion);

			Utils::writeLog("INFO. StartJVM -> AttachJVM starting...");
			AttachJVM();
			/*if(_env == 0)
			{
				writeLog("INFO. StartJVM -> AttachJVM starting...");
				AttachJVM();
			}
			else
			{
				writeLog("ERROR. StartJVM -> nVMs>0 but GetEnv failed.");
				CreateJVM();
			}*/
		}
		/// NO running instance
		else
		{
			Utils::writeLog("INFO. StartJVM -> nVMs=0");
			CreateJVM();
		}
	}
}

void JVM::StopJVM()
{
	///If attached, then  dettach from the VM
	if(GetThreadId(jvm) > 0)
	{
		Utils::writeLog("INFO. StopJVM -> DetachCurrentThread.");
		jvm->DetachCurrentThread();
	}
	///If created, then unload the Java VM (on program exit)
	else
	{
		Utils::writeLog("INFO. StopJVM -> DestroyJavaVM & FreeLibrary");
		jvm->DestroyJavaVM();
		FreeLibrary(hinstLib);
	}

	/// Destroy global references
	//env->DeleteGlobalRef(...);
}

void JVM::verifyEXEC()
{
	if (env->ExceptionCheck())
	{
		env->ExceptionDescribe(); 
		env->ExceptionClear();
	}
}

void JVM::test_InvokeMethod(JNIEnv* env)
{
	/// FIND class
	jclass myClass = (*env).FindClass("MyTest");
	
	/// CALL method
	if(myClass != nullptr)
	{
		std::cout << "Class found!";
		test_JavaCallback(env, myClass);
		test_CppCallback(env, myClass);
	}
	else
	{                                  
		std::cerr << "ERROR: class not found!";
	}
}

void JVM::test_JavaCallback(JNIEnv* env, jclass myClass)
{
	/// CALL method #1 (static function)
	jmethodID mid1 = env->GetStaticMethodID(myClass, "mymain", "()V");
	if(mid1 != nullptr)
	{
		env->CallStaticVoidMethod(myClass, mid1);                      
	}
	else
	{
		std::cerr << "ERROR: method not found !" << std::endl;
	}
		
	/// CALL method #2 (static function with integer)
	jmethodID mid2 = env->GetStaticMethodID(myClass, "mymain2", "(I)I");  
	if(mid2 == nullptr)
	{
		std::cerr << "ERROR: method it main2(int) not found !" << std::endl;
	}
	else
	{
		env->CallStaticVoidMethod(myClass, mid2, (jint)5);
	}

	/// CALL method #3 (static function with string array)

	jmethodID mid3 = env->GetStaticMethodID(myClass, "main", "([Ljava/lang/String;)V");
	if(mid3 == nullptr)
	{
		std::cerr << "ERROR: method not found !" << std::endl;
	}
	else
	{
		jobjectArray arr = env->NewObjectArray(
			5,											/// constructs java array of 5
			env->FindClass("java/lang/String"),			/// Strings
			env->NewStringUTF("str")					/// each initialized with value "str"
		);					
		env->SetObjectArrayElement(arr, 1, env->NewStringUTF("MYOWNSTRING"));	/// change an element
		env->CallStaticVoidMethod(myClass, mid3, arr);	/// call the method with the arr as argument.
		env->DeleteLocalRef(arr);						/// release the object
	}
}

void JVM::test_CppCallback(JNIEnv* env, jclass myClass)
{
	/// CREATE object and CALL method

	/// i. Register method 
	//JNINativeMethod methods[] { { "doTest", "()V", (void *)&test_CppMethod } };  /// mapping table
	auto method1 = &JVM::test_CppMethod;
	JNINativeMethod methods[] { { "doTest", "()V", (void *)&method1 } };  /// mapping table

	if(env->RegisterNatives(myClass, methods, 1) < 0)						/// register it
	{                        
		if(env->ExceptionOccurred())                                        /// verify if it's ok
			std::cerr << " OOOOOPS: exception when registreing naives" << std::endl;
		else
			std::cerr << " ERROR: problem when registreing naives" << std::endl;
	}

	/// ii. rerun object construction and method call
	jmethodID ctor = env->GetMethodID(myClass, "<init>", "()V");  /// FIND AN OBJECT CONSTRUCTOR 
	if(ctor == nullptr)
	{
		std::cerr << "ERROR: constructor not found !" << std::endl;
	}
	else
	{
		std::cout << "Object succesfully constructed !" << std::endl;
		jobject myo = env->NewObject(myClass, ctor);	          /// CREEATE OBJECT

		if(myo) {                                     /// IF OBJECT CREATED EXECUTE METHOD
			jmethodID show = env->GetMethodID(myClass, "showId", "()V");
			if(show == nullptr)
				std::cerr << "No showId method !!" << std::endl;
			else
				env->CallVoidMethod(myo, show);
		}
	}

	/// Callback from Java side
	/*MyTest
	{
		...
		public native void doTest();  // to be supplied in C++ trhough JNI

		public void showId() {
			System.out.println(uid);
			doTest();				// invoke the native method
		}
	}*/
}

void JVM::test_CppMethod(JNIEnv*e, jobject o) {
	std::cout << "C++callback activated" << std::endl;
	jfieldID f_uid = e->GetFieldID(e->GetObjectClass(o), "uid", "I");
	if (f_uid)
		 std::cout << "UID data member: " << e->GetIntField(o, f_uid) << std::endl;
	else 
		std::cout << "UID not found" << std::endl;
}