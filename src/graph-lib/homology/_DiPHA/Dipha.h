/*
* Pers.h
*
* Contains main() function for Persistent Homology
* Using Discrete Morse Theory
*/

#ifndef DIPHA_H
#define DIPHA_H
#pragma once

#include "include/dipha/includes.h"

class PH_DIPHA
{
public:
	/// Declaration of methods


	void PH_DIPHA::computePH(int a)
	{
		/// Initialize MPI (Message Passing Interface) standard for parallel execution
		MPI_Init(NULL, NULL);

		/// Define simplicial complex
		///dipha::inputs::abstract_weighted_cell_complex<int> wecellcomplex;
		dipha::inputs::weighted_explicit_complex wexcomplex;

		/// LOAD COMPLEX FROM FILE
		///	format: "weighted regular cell complex" in (co-)boundary matrix form (WEIGHTED_BOUNDARY_MATRIX):
		/// 0 if the file contains the boundaries of the cells, 1 if the file contains the coboundaries
		/// total number of cells n
		/// dimension of the cell complex d
		/// dimensions of the cells: dim(1) ... dim(n)
		/// floating point values of the cells: value(1) ... value(n)
		/// offsets of (co-)boundary entries of the cells: offset(1) ... offset(n). Each offset represents the first index of the (co-)boundary in the flattened (co-)boundary matrix.
		/// total number of non-zero entries m of the (co-)boundary matrix
		/// entries of the flattened (co-)boundary matrix: entry(1) ... entry(m)
		const std::string& input_filename = "D:\\Graph\\graph-suite\\graph-lib\\homology\\DiPHA\\test_data\\dual_explicit.complex";
		bool dualize = true;
		int64_t upper_dim = 1;
		const std::string& output_filename = "D:\\Graph\\graph-suite\\graph-lib\\homology\\DiPHA\\test_data\\output.bin";

		wexcomplex.load_binary(input_filename, upper_dim);

		if (dipha::globals::benchmark)
			dipha::mpi_utils::cout_if_root() << std::endl << "Number of cells in input: " << std::endl << wexcomplex.get_num_cells() << std::endl;
		dipha::data_structures::distributed_vector< int64_t > filtration_to_cell_map;
		dipha::data_structures::write_once_column_array reduced_columns;
		dipha::algorithms::compute_reduced_columns(wexcomplex, dualize, upper_dim, filtration_to_cell_map, reduced_columns);

		dipha::outputs::save_persistence_diagram(output_filename, wexcomplex, filtration_to_cell_map, reduced_columns, dualize, upper_dim);
	}

};

#endif