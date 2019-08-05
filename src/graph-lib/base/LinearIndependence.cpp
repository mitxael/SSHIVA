/// ***********************************************************
///	File	: Utils.cpp
///	About	: Gaussian Elimination's collection of variants
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
/// ******************************************************************

#include "LinearIndependence.h"		/// header file

/// Definition of methods

void swap_rows(std::vector<std::vector<double>> &A, std::size_t i, std::size_t k)
{
	// check indices
	assert(0 <= i);
	assert(i <= (A.size() - 1));

	assert(0 <= k);
	assert(k <= (A.size() - 1));

	for (std::size_t col = 0; col <= (A[0].size() - 1); ++col)
		std::swap(A[i][col], A[k][col]);
}

void divide_row(std::vector<std::vector<double>> &A, std::size_t i, double v)
{

	assert(0 <= i);
	assert(i <= (A.size() - 1));

	assert(v != 0);

	for (std::size_t col = 0; col <= (A[0].size() - 1); ++col)
		A[i][col] /= v;
}

/// in matrix A, add (row k*v) to (row i)
void add_multiple_row(std::vector<std::vector<double>>& A, std::size_t i, std::size_t k,double v)
{
	assert(0 <= i);
	assert(i <= (A.size() - 1));

	assert(0 <= k);
	assert(k <= (A.size() - 1));

	for (std::size_t col = 0; col <= (A[0].size() - 1); ++col)
		A[i][col] += v * A[k][col];
}

/// convert A to reduced row echelon form
void to_reduced_row_echelon_form(std::vector<std::vector<double>> &B)
{
	std::size_t lead = 0;

	for (std::size_t row = 0; row <= (B.size() - 1); ++row)
	{
		#ifdef _DEBUG
		if (verbosity > 2) {
			std::cout << "** Processing row " << row << ":" << std::endl;
			std::cout << "MCB before start: " << std::endl;
			for (std::size_t r = 0; r < B.size(); ++r) {
				for (std::size_t c = 0; c < B[r].size(); ++c)
					std::cout << B[r][c] << " ";
				std::cout << std::endl;
			}
		}
		#endif

		if (lead > (B[0].size() - 1))
			return;
		std::size_t i = row;
		while (B[i][lead] == 0)
		{
			++i;
			if (i > (B.size() - 1))
			{
				i = row;
				++lead;
				if (lead > (B[0].size() - 1))
					return;
			}
		}
		swap_rows(B, i, row);
		divide_row(B, row, B[row][lead]);
		for (i = 0; i <= (B.size() - 1); ++i)
		{
			if (i != row)
				add_multiple_row(B, i, row, -B[i][lead]);
		}
	}
}

bool gaussianElimination_dummy(std::vector<std::vector<int>> &A) {
	// This simple and dummy method to check that only one "1" is present on each column.
	int counter;			// keep track of "0"'s
	std::size_t ind = A.size() - 1; // test the last vector (against all previous ones)

	for (std::size_t col = 0; col < A[0].size(); col++) { // start looking for the first "1" in the last vector
		if (A[ind][col] == 1) { // if a "1" is found in the last vector...
			counter = 0;
			for (std::size_t row = ind - 1; (row >= 0 && row < 4294967295); row--) // check that no other "1" in found in the same column
				if (A[row][col] == 0)
					counter++;
			if (counter == ind)
			{
				#ifdef _DEBUG
				if (verbosity > 2) std::cout << "All vectors are linearly Independent" << std::endl;
				#endif

				return true;	// if one exclusive "1" was found at the end of the secondary loop, the last vector is linearly INDEPENDENT
			}

		}
	}
	
	#ifdef _DEBUG
	if (verbosity > 2) std::cout << "Two or more vectors are linearly Dependent" << std::endl;
	#endif

	return false;				// if no exclusive "1" was found at the end of the main loop, the last vector is linearly DEPENDENT
}

bool gaussianElimination_standard(std::vector<std::vector<int>> &A)
{
	double row_sum = 0;

	// Parse MCB from INT to DOUBLE
	std::vector<std::vector<double>> B;
	for (std::size_t i = 0; i < A.size(); i++) {
		std::vector<double> tmp(A[i].begin(), A[i].end());
		B.push_back(tmp);
		tmp.clear();
	}

	#ifdef _DEBUG
	if (verbosity > 2)
	{
		std::cout << "MCB before elimination: " << std::endl;
		for (std::size_t r = 0; r < B.size(); ++r) {
			for (std::size_t c = 0; c < B[r].size(); ++c)
				std::cout << B[r][c] << " ";
			std::cout << std::endl;
		}
	}
	#endif
	//perform gaussian elimination
	to_reduced_row_echelon_form(B);

	#ifdef _DEBUG
	if (verbosity > 2)
	{
		std::cout << "MCB after elimination:" << std::endl;
		for (std::size_t r = 0; r < B.size(); ++r) {
			for (std::size_t c = 0; c < B[r].size(); ++c)
				std::cout << B[r][c] << " ";
			std::cout << std::endl;
		}
	}
	#endif

	//check for "zero" rows
	for (std::size_t row = 0; row <= (B.size() - 1); ++row) {
		row_sum = 0;
		for (std::size_t col = 0; col < (B[0].size() - 1); col++) // adding all the values of the row
			row_sum += std::fabs(B[row][col]);					  // |absolute| as negative values can cause "false" zeros

		#ifdef _DEBUG
		if (verbosity > 2)
		{
			if (row_sum == 0) {
				std::cout << "Zero row!" << std::endl;
				std::cout << "--------------------------------------------------------" << std::endl;
				return false;	// if at least one "zero" row is found, cycles are NOT linearly-independent
			}
			std::cout << "--------------------------------------------------------" << std::endl;
		}
		#endif
	}

	return true;			// if no "zero" row is found, cycles are linearly-independent	
}

bool gaussianElimination_binary(std::vector<std::vector<int>> &A) {

	std::vector<std::vector<int>> B = A;
	std::size_t max_col = B[0].size() - 1;		//cols
	std::size_t max_row = B.size() -1 ;		//rows
	std::size_t col_pivot = 0;

	//perform "gaussian" elimination
	for (std::size_t row = 0; row <= max_row; ++row) //for every column
	{
		#ifdef _DEBUG
		if (verbosity > 2) {
			std::cout << "** Processing row " << row << ":" << std::endl;
			std::cout << "MCB before start: " << std::endl;
			for (std::size_t r = 0; r < B.size(); ++r) {
				for (std::size_t c = 0; c < B[r].size(); ++c)
					std::cout << B[r][c] << " ";
				std::cout << std::endl;
			}
		}
		#endif

		if (col_pivot > max_col)
			break;
				
		// Search for first row having "1" in column "lead"
		std::size_t i = row; //Start searching from row "i"=row
		while (B[i][col_pivot] == 0)
		{
			++i;
			if (i > max_row) // no "1" was found across rows
			{
				i = row;		// set "i" back to "row"
				++col_pivot;	// move to next column
				if (col_pivot > max_col) //if no more columns
					break;
			}
		}
		
		// swap "i" and "row" rows
		if ((0 <= i) && (i <= max_row) && (0 <= row) && (row <= max_row)) {
			#ifdef _DEBUG
			if (verbosity > 2) std::cout << "Swap row " << i << " with " << row << "." << std::endl;
			#endif
			for (std::size_t col = 0; col <= max_col; ++col) {
				std::swap(B[i][col], B[row][col]);
			}
		}
		
		// perform XOR in "k-th" row against "row-th" row
		for (int k = 0; k <= max_row; ++k) //k=col_pivot?
		{
			if ( (k != row)   ) // && (B[k][row] != 0) //Apply only to "1"s
			{
				int check_zero_row = 0;
				#ifdef _DEBUG
				if (verbosity > 2) std::cout << "XOR row " << k << " with " << row << "." << std::endl;
				#endif
				for (std::size_t j = 0; j <= max_col; j++) {
					B[k][j] = B[k][j] ^ B[row][j];				//XOR: Use "1" to set other "1" rows to 0 (and 0 to 1)
					check_zero_row = check_zero_row | B[k][j];	//OR: Result will be "0" only if the entire row is "0"
				}
				if (check_zero_row == 0 )
				{
					#ifdef _DEBUG
					if (verbosity > 2) std::cout << "Zero row found at row" << k << "! " << std::endl;
					#endif
					return false;
				}
			}
		}
	}

	#ifdef _DEBUG
	if (verbosity > 2) {
		std::cout << "MCB after elimination: " << std::endl;
		for (std::size_t r = 0; r <= max_row; ++r) {
			for (std::size_t c = 0; c < B[r].size(); ++c)
				std::cout << B[r][c] << " ";
			std::cout << std::endl;
		}
	}
	#endif

	//check if LastRow is "zero"
	/*int row_sum = 0;
	for (std::size_t col = 0; col <= max_col; col++)									// adding all the values of the row
		row_sum += std::fabs(B[max_row][col]);											// |absolute| as negative values can cause "false" zeros

	if (row_sum == 0) {
		if (verbosity > 2) {
			std::cout << "Zero row!" << std::endl;
			std::cout << "--------------------------------------------------------" << std::endl;
		}
		return false;																// if at least one "zero" row is found, cycles are NOT linearly-independent
	}*/

	#ifdef _DEBUG
	if (verbosity > 2) std::cout << "--------------------------------------------------------" << std::endl;
	#endif

	return true;																		// if no "zero" row is found, cycles are linearly-independent
}

bool gaussianElimination_bitset_mod2(std::vector<std::vector<int>> &A)
{
	//BITSET VERSION
	const int N = 8;
	int n = static_cast<int>(A.size());
	int m = static_cast<int>(A[0].size());

	//CONVERT A to a
	std::vector <std::bitset<N>> a(n);
	std::string temp_str;
	for (int i = 0; i < n; ++i)
	{
		for (int j = m; j >= 0; --j)
		{
			temp_str = temp_str + std::to_string(A[i][j]);
		}
		std::bitset<N> temp_bin(temp_str);
		a[i] = temp_bin;
		temp_str.clear();
	}

	//SOLVE
	for (int col = 0, row = 0; col < m && row < n; ++col) {
		for (int i = row; i < n; ++i)
			if (a[i][col]) {
				swap(a[i], a[row]);
				break;
			}
		if (!a[row][col])
			continue;
		for (int i = 0; i < n; ++i)
			if (i != row && a[i][col])
				a[i] ^= a[row];
		++row;
	}

	return (a[n-1]==0) ? false : true;
}

bool gaussianElimination_mod2(const std::vector<std::vector<ut>> &A)
{
	int n = static_cast<int>(A.size());
	int m = static_cast<int>(A[0].size());
	std::vector <std::vector<ut>> a = A;

	/// CONVERT to boolean, if required by parameter
	/*bool convertToBool = false;
	if (convertToBool)
	{
		std::vector <std::vector<bool>> a(n);
		std::string temp_str;
		for (int i = 0; i < n; ++i)
			for (int j = 0; j < m; ++j)
				a[i].push_back(A[i][j]);
	}*/

	/// Work on all columns and all but last row
	for (int col = 0, row = 0; col < m && row < n -1 ; ++col)
	{
		/// Look for column holding the first "1" in the "pivot" row
		for (int i = row; i < n; ++i)
		{
			/// If a "1" is found in the row, swap it with the first inspected row
			if (a[i][col] == 1) {
				std::swap(a[i], a[row]);
				break;
			}
		}

		/// If the "pivot" row doesn't have any "1", jump to next column
		if (a[row][col] == 0)
			continue;

		/// XOR the pivot "row" against all other rows
		for (int i = row; i < n; ++i)
		{
			/// XOR only against other rows with a "1" in the pivot column (to convert the "1" into "0", and so the entire column)
			if (i != row && a[i][col] == 1)
			{
				for (int j = row; j < m; ++j)
				{
					a[i][j] = a[i][j] ^ a[row][j];
				}
			}
		}
		++row;	/// row increases only when a "1" column was found and processed
	}

	int zerorow = 0;
	for (int k = 0; k < m; ++k)
		zerorow |= a[n-1][k];	/// cumulative OR

	#ifdef _DEBUG
	if (verbosity > 2 && zerorow == 0)
		std::cout << "*** Cycle is linearly independent mod 2." << std::endl;
	#endif

	return zerorow;
}