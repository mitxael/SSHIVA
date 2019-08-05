/*
* Perseus.cpp
*
* Contains main() function for Persistent Homology
* Using Discrete Morse Theory
*/


#include "Perseus.h"

#include "Cells/Cell.h"
#include "Cells/Cube.h"
#include "Cells/PCell.h"
#include "Cells/Point.h"
#include "Cells/Simplex.h"
#include "Complexes/Complex.h"
#include "Complexes/CToplex.h"
#include "Complexes/DenseCToplex.h"
#include "Complexes/Morse.h"
#include "Complexes/PComplex.h"
#include "Complexes/RIPS.h"
#include "Complexes/SToplex.h"

typedef int BC;
typedef int CC;
void PH_PERSEUS::computePH(int a)
{
	/// Create Complex
	MComplex<CC, BC>* ccomp;
	std::string infile = "path\\to\\input";
	std::string outfile = "path\\to\\output";
	std::ifstream inf(infile.c_str(), std::ifstream::in);
	SToplex<CC, double, BC> stop = SToplex<CC, double, BC>();
	bool hascap = false;
	bool hasbirths = true;
	int complexType = 0;	///0:NMST; 1:ST
	if (inf.good())
	{
		std::pair<num, bool> iores;
		/*if (complexType == 0)
		{
			/// BUILD cell complex for NMST non-manifold simplicial toplex from file. the input TEXT file inf must have
			/// the following format:
			/// <dimension of embedding space n>
			/// <dim of simplex, m> <x_11,...,x_1n> <x_21,...,x_2n> ... <x_(m+1)1,...x_(m+1)n> <birth time>
			/// : [ similarly for other simplices, each with its own dim m ]
			iores = stop.makeFromFile(inf, hasbirths, false, hascap);
		}
		else
		{
			/// BUILD cell complex from simplicial toplex information. the TEXT file inf
			/// should have the following format:
			/// <simplex dimension m>
			/// <ambient space dimension n>
			/// <x_11,...,x_1n> <x_21,...,x_2n> ... <x_(m+1)1,...x_(m+1)n> <birth time> { i.e. vertex coords and birth time }
			/// : [same as above for next simplex]
			/// : [ etc ]
			iores = stop.makeFromFile(inf, hasbirths, true, hascap);
		}

		/// COMPUTE HOMOLOGY
		if (iores.second != false)
		{
			if (FLOWTALK) std::cout << "\nRead " << iores.first << " top simplices from Input File ";
			if (MAKEBPS) std::cin.get();

			if (FLOWTALK) std::cout << "\nWriting Cell Complex From Non-manifold simplicial Complex";
			if (MAKEBPS) std::cin.get();

			//cout<<stop;
			stop.writeComplex(*ccomp);

			//ccomp.checkComplex();
			if (FLOWTALK) std::cout << "\nDone!";
			if (MAKEBPS) std::cin.get();
		}


		if (FLOWTALK)
		{
			std::cout << " Complex stored with " << ccomp->size() << " cells!";
			if (MAKEBPS) std::cin.get();
		}

		/// optimize!!
		ccomp->hyperq = true;
		std::string eng = "r";
		bool savegens = true;
		if (eng == "r")
		{
			ccomp->ReduceAndUpdate(ccomp, savegens, 0.2, false, std::cout);
		}
		else if (eng == "c")
		{
			ccomp->CoreduceAndUpdate(ccomp, savegens, 0.2, false, std::cout);
		}
		else if (eng == "s") // skip (co)reductions altogether, just compute persistence
		{

		}
		else
		{
			ccomp->AlternateAndUpdate(ccomp, savegens, 0.2, false, true, std::cout);
		}

		PComplex<CC, BC> pcomp;
		if (FLOWTALK) std::cout << "\nComputing Persistence Intervals!";
		if (MAKEBPS) { std::cout << "... ";    std::cin.get(); }
		bool truncate = false;
		pcomp.COMPUTE_INTERVALS(*ccomp, savegens, truncate);

		//ccomp.showInts();
		pcomp.makeOutputFiles(outfile);
		//pcomp.showBetti();

		delete ccomp;
		//delete ccomp;

		if (FLOWTALK) std::cout << "\n\nDone!!! Please consult [" << outfile << "*.txt] for results.\n\n";*/
	}

	inf.close();
}
