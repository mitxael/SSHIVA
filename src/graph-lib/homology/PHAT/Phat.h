/*
* PHAT.h
*
* Contains main() function for Persistent Homology
* Using Discrete Morse Theory
*/

#ifndef PHAT_H
#define PHAT_H
#pragma once

#include "include/phat/compute_persistence_pairs.h"			/// PHAT
#include "include/phat/boundary_matrix.h"					/// PHAT
#include "include/phat/persistence_pairs.h"					/// PHAT
#include "include/phat/algorithms/standard_reduction.h"		/// PHAT


class PH_PHAT
{
public:
	/// Declaration of methods

	void PH_PHAT::computePH(Graph& G, int a)
	{
		/// boundary_matrix - ascii:
		/// The file represents the filtration of the cell complex, containing one cell per line (empty lines and lines starting with "#" are ignored).
		/// A cell is given by a sequence of integers, separated by spaces, where the first integer denotes the dimension of the cell, 
		/// and all following integers give the indices of the cells that form its boundary (the index of a cell is its position in the filtration, 
		/// starting with 0). A sample file single_triangle.dat can be found in the examples folder.

		/*
		///std::cout << "We will build an ordered boundary matrix of this simplicial complex consisting of a single triangle: " << std::endl;
		///std::cout << std::endl;
		///std::cout << " 3" << std::endl;
		///std::cout << " |\\" << std::endl;
		///std::cout << " | \\" << std::endl;
		///std::cout << " |  \\" << std::endl;
		///std::cout << " |   \\ 4" << std::endl;
		///std::cout << "5|    \\" << std::endl;
		///std::cout << " |     \\" << std::endl;
		///std::cout << " |  6   \\" << std::endl;
		///std::cout << " |       \\" << std::endl;
		///std::cout << " |________\\" << std::endl;
		///std::cout << " 0    2    1" << std::endl;

		/// first define a boundary matrix with the chosen internal representation
		phat::boundary_matrix< phat::vector_vector > boundary_matrix;

		/// set the number of columns (has to be 7 since we have 7 simplices)
		boundary_matrix.set_num_cols(7);

		/// set the dimension of the cell that a column represents:
		boundary_matrix.set_dim(0, 0);
		boundary_matrix.set_dim(1, 0);
		boundary_matrix.set_dim(2, 1);
		boundary_matrix.set_dim(3, 0);
		boundary_matrix.set_dim(4, 1);
		boundary_matrix.set_dim(5, 1);
		boundary_matrix.set_dim(6, 2);

		/// set the respective columns -- the columns entries have to be sorted
		std::vector< phat::index > temp_col;

		boundary_matrix.set_col(0, temp_col);	///0
		boundary_matrix.set_col(1, temp_col);	///1
		temp_col.push_back(0);
		temp_col.push_back(1);
		boundary_matrix.set_col(2, temp_col);	///0 1
		temp_col.clear();

		boundary_matrix.set_col(3, temp_col);	///3
		temp_col.push_back(1);
		temp_col.push_back(3);
		boundary_matrix.set_col(4, temp_col);	///1 3
		temp_col.clear();

		temp_col.push_back(0);
		temp_col.push_back(3);
		boundary_matrix.set_col(5, temp_col);	///0 3
		temp_col.clear();

		temp_col.push_back(2);
		temp_col.push_back(4);
		temp_col.push_back(5);
		boundary_matrix.set_col(6, temp_col);	///2 4 5
		temp_col.clear();*/

		///TODO: All added with a different filtration value => 24!!
		///ALL CLIQUES COMPLEX
		phat::boundary_matrix< phat::vector_vector > boundary_matrix;
		buildAllCliquesComplex(G, boundary_matrix);

		/// print some information of the boundary matrix:
		/*std::cout << std::endl;
		std::cout << "The boundary matrix has " << boundary_matrix.get_num_cols() << " columns: " << std::endl;
		for (phat::index col_idx = 0; col_idx < boundary_matrix.get_num_cols(); col_idx++) {
			std::cout << "Column " << col_idx << " represents a cell of dimension " << (int)boundary_matrix.get_dim(col_idx) << ". ";
			if (!boundary_matrix.is_empty(col_idx)) {
				std::vector< phat::index > temp_col;
				boundary_matrix.get_col(col_idx, temp_col);
				std::cout << "Its boundary consists of the cells";
				for (phat::index idx = 0; idx < (phat::index)temp_col.size(); idx++)
					std::cout << " " << temp_col[idx];
			}
			std::cout << std::endl;
		}
		std::cout << "Overall, the boundary matrix has " << boundary_matrix.get_num_entries() << " entries." << std::endl;*/


		/// define the object to hold the resulting persistence pairs
		phat::persistence_pairs pairs;

		/// choose an algorithm (choice affects performance) and compute the persistence pair
		/// (modifies boundary_matrix)
		phat::compute_persistence_pairs< phat::twist_reduction >(pairs, boundary_matrix);

		/// sort the persistence pairs by birth index 
		pairs.sort();

		/// print the pairs:
		std::cout << "Intervals: " << std::endl;
		for (phat::index idx = 0; idx < pairs.get_num_pairs(); idx++)
		{
			std::cout << "\t[" << pairs.get_pair(idx).first << ", " << pairs.get_pair(idx).second << ")" << std::endl;

		}
		
		/// Add barcodes to graph
		barcodesIntoGraph(G, pairs);
	}

	void PH_PHAT::barcodesIntoGraph(Graph& G, phat::persistence_pairs pairs)
	{
		
	}

	void PH_PHAT::compute_pairing(bool use_binary, bool verbose, bool dualize)
	{
		/// Files
		std::string file_output = "D:\\Graph\\graph-suite\\graph-lib\\homology\\PHAT\\_examples\\mitxael.dat";	///persistance pairs
		std::string file_boundarymatrix = "D:\\Graph\\graph-suite\\graph-lib\\homology\\PHAT\\_examples\\single_triangle.dat";
		std::string file_pairs = "D:\\Graph\\graph-suite\\graph-lib\\homology\\PHAT\\_examples\\single_triangle_persistence_pairs.dat";

		/// Compute persistance pairs
		phat::boundary_matrix<phat::vector_vector> boundary_matrix;
		bool read_successful;

		double read_timer = omp_get_wtime();
		if (use_binary) {
			std::cout << "Reading input file " << file_boundarymatrix << " in binary mode" << std::endl;
			read_successful = boundary_matrix.load_binary(file_boundarymatrix);
		}
		else {
			std::cout << "Reading input file " << file_boundarymatrix << " in ascii mode" << std::endl;
			read_successful = boundary_matrix.load_ascii(file_boundarymatrix);
		}
		double read_time = omp_get_wtime() - read_timer;
		double read_time_rounded = floor(read_time * 10.0 + 0.5) / 10.0;
		std::cout << "Reading input file took " << std::setiosflags(std::ios::fixed) << std::setiosflags(std::ios::showpoint) << std::setprecision(1) << read_time_rounded << "s" << std::endl;

		if (!read_successful) {
			std::cerr << "Error opening file " << file_boundarymatrix << std::endl;
		}

		phat::index num_cols = boundary_matrix.get_num_cols();

		if (dualize) {
			double dualize_timer = omp_get_wtime();
			std::cout << "Dualizing ..." << std::endl;
			phat::dualize(boundary_matrix);
			double dualize_time = omp_get_wtime() - dualize_timer;
			double dualize_time_rounded = floor(dualize_time * 10.0 + 0.5) / 10.0;
			std::cout << "Dualizing took " << std::setiosflags(std::ios::fixed) << std::setiosflags(std::ios::showpoint) << std::setprecision(1) << dualize_time_rounded << "s" << std::endl;
		}

		double pairs_timer = omp_get_wtime();
		phat::persistence_pairs pairs;
		std::cout << "Computing persistence pairs ..." << std::endl;
		phat::compute_persistence_pairs <phat::standard_reduction>(pairs, boundary_matrix);
		double pairs_time = omp_get_wtime() - pairs_timer;
		double pairs_time_rounded = floor(pairs_time * 10.0 + 0.5) / 10.0;
		std::cout << "Computing persistence pairs took " << std::setiosflags(std::ios::fixed) << std::setiosflags(std::ios::showpoint) << std::setprecision(1) << pairs_time_rounded << "s" << std::endl;

		if (dualize) dualize_persistence_pairs(pairs, num_cols);


		double write_timer = omp_get_wtime();
		if (use_binary) {
			std::cout << "Writing output file " << file_output << " in binary mode ..." << std::endl;
			pairs.save_binary(file_output);
		}
		else {
			std::cout << "Writing output file " << file_output << " in ascii mode ..." << std::endl;
			pairs.save_ascii(file_output);
		}
		double write_time = omp_get_wtime() - write_timer;
		double write_time_rounded = floor(write_time * 10.0 + 0.5) / 10.0;
		std::cout << "Writing output file took " << std::setiosflags(std::ios::fixed) << std::setiosflags(std::ios::showpoint) << std::setprecision(1) << write_time_rounded << "s" << std::endl;
	}

	void PH_PHAT::buildAllCliquesComplex(Graph& G, phat::boundary_matrix< phat::vector_vector >& boundary_matrix)
	{
		/// FLAT PH_filtration: std::vector<std::pair<id_simplex, id_boundary[]>>
		std::vector<std::pair<int, std::vector<int>>> flatFiltration;

		/// Determine number of simplices
		int n_simplices = 0;
		for (auto it = G.PH_Filtration.begin(); it != G.PH_Filtration.end(); ++it)
		{
			for (auto cliqueP : it->second)
			{
				n_simplices++;
			}
		}
		boundary_matrix.set_num_cols(n_simplices);
		
		/// For all cliques (in filtration order)
		std::vector<phat::index> col_simplex;
		int id_simplex = 0;
		for (auto it = G.PH_Filtration.begin(); it != G.PH_Filtration.end(); ++it)
		{
			const int &rank = it->first;
			const auto &rankCliques = (it->second);

			/// Add cliques with size of given "rank"
			for (auto cliqueP : rankCliques)
			{
				std::vector<int> &clique = cliqueP.first;
				double weight = cliqueP.second;
				int dim_simplex = clique.size()-1;
				//int k = (int)clique.size();

				/// Add simplex and dimension
				boundary_matrix.set_dim(id_simplex, dim_simplex);

				/// Add boundary of [dim]-simplex (i.e. signed sum of its [dim-1]-simplices)
				/// should add the id_simplex of [dim-1]-simplices which form it
				if (dim_simplex > 0)
				{
					col_simplex.clear();
					for (auto vertex : clique)
					{
						/// Create new clique as the copy of source clique BUT excluding the current "vertex"
						std::vector<int> boundary;
						std::copy_if(clique.begin(), clique.end(), std::back_inserter(boundary), [vertex](const int& elem) { return elem != vertex; });	

						/// find the id_simplex of "boundary"
						int id_boundary = 0;
						int found = 0;
						for (auto it = flatFiltration.rbegin(); it != flatFiltration.rend(); ++it)
						{
							auto lookup_id = (*it).first;
							auto lookup_value = (*it).second;
							
							if(lookup_value.size() == dim_simplex)
							{
								///if(boundary == lookup_value)
								if( compareVectors(boundary,lookup_value) == true )
								{
									col_simplex.push_back(lookup_id);
									found++;
								}
							}
							if(found == dim_simplex)
								break;
						}	
					}
					std::sort (col_simplex.begin(), col_simplex.end());
					boundary_matrix.set_col(id_simplex, col_simplex);
				}

				/// Add simplex to flatFiltration
				flatFiltration.push_back(std::make_pair(id_simplex, clique));

				///increase id
				id_simplex++;
			}
		}
	}

	/*void PH_PHAT::buildAllCliquesComplex(Graph& G, phat::boundary_matrix< phat::vector_vector >& boundary_matrix)
	{
		/// FLAT PH_filtration: std::vector<std::pair<id_simplex, id_boundary[]>>
		std::vector<std::pair<int, std::vector<int>>> flatFiltration;

		/// Determine number of simplices
		int n_simplices = 0;
		for (auto it = G.PH_Filtration.begin(); it != G.PH_Filtration.end(); ++it)
		{
			for (auto cliqueP : it->second)
			{
				/// Add simplex to flatFiltration
				int id_simplex = n_simplices;
				std::vector<int> simplex = cliqueP.first;
				flatFiltration.push_back(std::make_pair(id_simplex, simplex));

				n_simplices++;
			}
		}
		boundary_matrix.set_num_cols(n_simplices);

		/// Add simplices and dimension
		n_simplices = 0;
		for (auto it = G.PH_Filtration.begin(); it != G.PH_Filtration.end(); ++it)
		{
			for (auto cliqueP : it->second)
			{
				int id_simplex = n_simplices;
				int dim_simplex = cliqueP.first.size()-1;
				boundary_matrix.set_dim(id_simplex, dim_simplex);

				n_simplices++;
			}
		}
		
		/// For all cliques (in filtration order)
		std::vector<phat::index> col_simplex;
		int id_simplex = 0;
		for (auto it = G.PH_Filtration.begin(); it != G.PH_Filtration.end(); ++it)
		{
			const int &rank = it->first;
			const auto &rankCliques = (it->second);

			/// Add cliques with size of given "rank"
			for (auto cliqueP : rankCliques)
			{
				std::vector<int> &clique = cliqueP.first;
				double weight = cliqueP.second;
				int dim_simplex = clique.size()-1;
				//int k = (int)clique.size();			

				/// Add boundary of [dim]-simplex (i.e. signed sum of its [dim-1]-simplices)
				/// should add the id_simplex of [dim-1]-simplices which form it
				if (dim_simplex > 0)
				{
					col_simplex.clear();
					for (auto vertex : clique)
					{
						/// Create new clique as the copy of source clique BUT excluding the current "vertex"
						std::vector<int> boundary;
						std::copy_if(clique.begin(), clique.end(), std::back_inserter(boundary), [vertex](const int& elem) { return elem != vertex; });	

						/// find the id_simplex of "boundary"
						int id_boundary = 0;
						int found = 0;
						for (auto it = flatFiltration.rbegin(); it != flatFiltration.rend(); ++it)
						{
							auto lookup_id = (*it).first;
							auto lookup_value = (*it).second;
							
							if(lookup_value.size() == dim_simplex)
							{
								///if(boundary == lookup_value)
								if( compareVectors(boundary,lookup_value) == true )
								{
									col_simplex.push_back(lookup_id);
									found++;
								}
							}
							if(found == dim_simplex)
								break;
						}	
					}
					std::sort (col_simplex.begin(), col_simplex.end());
					boundary_matrix.set_col(id_simplex, col_simplex);
				}

				///increase id
				id_simplex++;
			}
		}
		bool end = true;
	}*/

	bool PH_PHAT::compareVectors(std::vector<int> vec1, std::vector<int> vec2)
	{
		bool elemFound;
		if(vec1.size() == vec2.size())
		{
			for (auto elem1 : vec1)
			{
				elemFound = false;
				for (auto elem2 : vec2)
				{
					if(elem1 == elem2)
					{
						elemFound = true;
						break;
					}
				}
				if(elemFound == false)
					return false;
			}
			return true;
		}
		else
			return false;
	}
};

#endif