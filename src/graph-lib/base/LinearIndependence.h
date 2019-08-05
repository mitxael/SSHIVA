/// ***********************************************************
///	File	: Utils.cpp
///	About	: Header of Linear Independence class
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
/// ******************************************************************

#ifndef LINEARINDEPENDENCE_H
#define LINEARINDEPENDENCE_H
#pragma once

#include "Graph.h"					/// Include the "extern" variables defined there

#include <vector>					/// std::vector
#include <cassert>					/// assert expressions
#include <iostream>					/// std::cin std::cout
#include <functional>				/// 
#include <bitset>					/// bitwise support

	// Declaration of member functions
	void swap_rows(std::vector<std::vector<double>> &A, std::size_t i, std::size_t k);					/// Swap rows i and k of a matrix A
	void divide_row(std::vector<std::vector<double>> &A, std::size_t i, double v);						/// divide row i of matrix A by v
	void add_multiple_row(std::vector<std::vector<double>>& A, std::size_t i, std::size_t k, double v); /// in matrix A, add (row k*v) to (row i)
	void to_reduced_row_echelon_form(std::vector<std::vector<double>> &B);								/// convert A to reduced row echelon form
	bool gaussianElimination_dummy(std::vector<std::vector<int>> &A);		/// Check if the LAST row of the matrix is "linearly independent" with respect to other rows
	bool gaussianElimination_standard(std::vector<std::vector<int>> &A);	/// Check if a matrix is a basis of "linearly independent" vectors
	bool gaussianElimination_binary(std::vector<std::vector<int>> &A);		/// Check if the LAST row of the matrix is "linearly independent" with respect to other rows
	bool gaussianElimination_bitset_mod2(std::vector<std::vector<int>> &A);	/// Check if a set of incidence "std::vector<std::bitset<N>>" is linearly independent
	bool gaussianElimination_mod2(const std::vector<std::vector<ut>> &A);		/// Check if a set of incidence vectors is linearly independent

#endif