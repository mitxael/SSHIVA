/// ******************************************************************
///	File	: Javaplex.h
///	About	: Header of Javaplex class
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// -  
/// ******************************************************************

#ifndef JAVAPLEX_H
#define JAVAPLEX_H
#pragma once
#pragma comment(lib,"advapi32.lib")									/// Required for JVM

#include "../../base/Graph.h"											/// Graph's header file
#include "../../base/Utils.h"											/// Utils's header file
#include "../../cliques/BronKerbosch.h"									/// Bron-Kerbosch's header file
#include "../../cycles/Horton.h"										/// Horton's header file

#include <jni.h>													/// Java Native Interface
#include <iostream>													/// Standard input/output stream objects
#include <map>														/// Map and Multimap associative containers
#include <unordered_map>											/// Unordered map and multimap associative containers
#include <list>		
#include <set>	
#include <unordered_set>
//#include <thread>													/// Not supported wehn compiling with CLR
//#include <ppl.h>													/// MS Parallel Patterns Library (concurrency::parallel_for) (not supported with /CRL)
#include <omp.h>													/// Open Multi Processing

/// Declaration of class, and its member variables and functions
	class PH_Javaplex
	{
	public:
		/// Declaration of methods
		
		jobject j_buildComplex(Graph &G, int type, int dimension);
		void j_buildGraphComplex(Graph &G, jobject &stream);					/// Build Simplicial Complex from vertices and edges
		void j_buildMinCyclesComplex(Graph &G, jobject &stream);				/// Build Simplicial Complex from Minimum Cycles Basis
		void j_buildMaxCliqueComplex(Graph &G, jobject &stream);				/// Build Simplicial Complex from Maximal Cliques Set
		void j_buildAllCliquesComplex(Graph &G, jobject &stream, bool autoAlgorithm, int dimension);	/// Build Simplicial Complex from All Cliques
		std::pair<jobject, jobject> computePH(Graph &G, jobject &stream);
		std::vector<int> getBettiNumbers(std::pair<jobject, jobject> &PH, bool onlyInfinite, bool getSequence);
		void getIntervalsAndGenerators(Graph &G, std::pair<jobject, jobject> &PH, int dimension, int displayMode, bool onlyInfinite, bool saveCycles);
		void holesIntoGraph(Graph &G, jobject annotatedBarcodes);
		void barcodesIntoGraph(Graph &G, jobject annotatedBarcodes);

		jobject createExplicitSimplexStream();								/// Creates a new Filtered Chain Complex
		void addVertex(jobject stream, const int &vertex, double filtrationValue=0);	/// Add a vertex to the complex (with either zero or given filtration)
		void addElement(jobject stream, const std::vector<int> &basisElements, double filtrationValue=0);	///Add a simplex with the provided vertices (with either zero or given filtration)
		void finalizeStream(jobject stream);								/// Prepare the stream for use by a consumer (e.g. PersistentHomology)
		int getSize(jobject stream);

		jobject getModularSimplicialAlgorithm(int maxDimension, int primeCoefficient);	/// Return a Simplicial Persistence Algorithm over the finite field "Z/pZ"
		jobject computeIntervals(jobject persistence, jobject stream);		/// Calculate the Persistent Homology Intervals of a given stream
		jobject getInfiniteIntervals(jobject barcodes);
		std::string getBettiNumbers(jobject barcodes);						/// Return the set of Betti barcodes = Ranks of homology groups = Number of holes in each dimension
		std::vector<int> getBettiSequence(jobject barcodes);				/// Return an array indicating the cardinality of the Betti barcodes at each dimension
		void ensureAllFaces(jobject stream);								/// Adds all faces missing to form the simplicial complex
		bool removeElementIfPresent(jobject stream, const std::vector<int> &basisElements);	/// fnc(0:(dimension + 1))	removes the k-simplex to respect k-1
		bool validateVerbose(jobject stream);								/// Check if the stream is a simplicial complex
		std::map<std::vector<std::vector<int>>, double> createBoundaryMatrixAsDoubleSum(jobject stream);

		jobject computeAnnotatedIntervals(jobject persistence, jobject stream);	/// Computes the AUGMENTED persistence intervals for a given stream (includes Holes' elements)
		std::string getAnnotatedIntervals_AsString(jobject annotatedBarcodes);
		std::map<int, int> getBettiNumbersMap(jobject annotatedBarcodes);	/// **TO BE REMOVED**
		std::vector<int> getDimensions(jobject annotatedBarcodes);			/// Return the set of dimensions at which there are intervals.
		
		std::multimap<std::string, std::string> getIntervalGeneratorPairsAtDimension(jobject annotatedBarcodes, int dim);	/// Return the set of interval-generator pairs at a specified dimension
		std::string getIntervalGeneratorPairsAtDimension_AsString(jobject annotatedBarcodes, int dim);	/// Return the set of interval-generator pairs at a specified dimension
		jobject getIntervalsAtDimension(jobject annotatedBarcodes, int dim);
		std::vector<std::vector<std::vector<int>>> getGeneratorsAtDimension(jobject annotatedBarcodes, int dim);	/// Return the set of generators at the given dimension
		std::vector<std::pair<double, double>> getEndPoints(jobject annotatedBarcodes, int dim, bool skipInfIntervals); /// Return an array containing the set of endpoint of the intervals

		std::string convert_jstring2string(jstring source);
		std::vector<std::vector<int>> convert_generator2vector(std::string generatorS);
		jintArray convert_intVector2jintArray(const std::vector<int> &basisElements);
		std::vector<int> convert_jintArray2intVector(const jintArray &basisElements);
		jdoubleArray convert_doubleVector2jdoubleArray(const std::vector<double> &basisElements);
		std::vector<double> convert_jdoubleArray2doubleVector(const jdoubleArray &basisElements);

		jobject Plex4();													/// Javaplex's API wrapper
		jobject Plex4_createExplicitSimplexStream(jobject plex4);
		jobject Plex4_getModularSimplicialAlgorithm(jobject plex4);

	};

	/// Declaration of class, and its member variables and functions
	class JVM
	{
	public:
		/// Declaration of variables
		/*HINSTANCE hinstLib;
		HMODULE	jvmHandle;
		JavaVM *jvm(0);								/// Pointer to the JVM (Java Virtual Machine)
		JNIEnv *env(0);								/// Pointer to native method interface
		jint JavaVersion = JNI_VERSION_1_6;			/// Minimum Java version required*/

		/// Declaration of methods
		JNIEnv* CreateJVM();												/// Create JVM to execute Java code
		JNIEnv* AttachJVM();												/// Attach JVM to execute Java code
		void StartJVM();													/// Determine whether to Create or Attach the JVM
		void StopJVM();														/// Destroy/Detach the JVM
		void verifyEXEC();													/// Check for exceptions on the JVM environment
		jsize GetnVMs();
		void test_InvokeMethod(JNIEnv* env);								/// Test Invoke() method
		void test_JavaCallback(JNIEnv* env, jclass myClass);				/// Test callback FROM Java
		void test_CppCallback(JNIEnv* env, jclass myClass);					/// Test callback TO Java
		void test_CppMethod(JNIEnv*e, jobject o);							/// Dummy method for test_CppCallback()
	};

#endif