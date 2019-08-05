///
/// Gudhi.h
///

#ifndef GUDHI_H
#define GUDHI_H
#pragma once


#include "include/gudhi/Simplex_tree.h"
#include "include/gudhi/Persistent_cohomology.h"
#include "include/gudhi/Persistence_intervals.h"

#include "../../base/Graph.h"											/// Graph's header file
#include "../../base/Utils.h"											/// Utils's header file
#include "../../cliques/BronKerbosch.h"									/// Bron-Kerbosch's header file
#include "../../cycles/Horton.h"										/// Horton's header file

using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Gudhi::Simplex_tree<>, Gudhi::persistent_cohomology::Field_Zp >;
using Persistence_intervals = Gudhi::Persistence_representations::Persistence_intervals;
using Simplex = std::vector< Gudhi::Simplex_tree<>::Vertex_handle >;
using Barcodes = std::vector<Persistent_cohomology::Persistent_interval>;

	class PH_GUDHI
	{
	public:
		static void computePH(Graph& G, bool onlyInfinite, int dimension, bool saveCycles);
		static void complexFromFiltration(Graph& G,  Gudhi::Simplex_tree<>& complex);
		static void buildGraphComplex(Graph& G, Gudhi::Simplex_tree<>& st);
		static void buildCliqueComplex(Graph& G,  Gudhi::Simplex_tree<>& complex);
		static int intervalsFromFile(int a);
		static void barcodesIntoGraph(Graph &G, Persistent_cohomology& pcoh);
		static void holesIntoGraph(Graph &G, Persistent_cohomology& ph);
		static double bottleneckDistance(std::string file1, std::string file2, double tolerance);
	};


#endif