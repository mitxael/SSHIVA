/*
 * PComplex.h
 *
 *  Persistent Complex
 *      Author: Vidit
 */

#ifndef PCOMPLEX_H_
#define PCOMPLEX_H_

#include "../Global/Global.h"
#include "../Cells/PCell.h"
#include "Complex.h"
#include "Morse.h"

// vector of persistent cells!
#define PCVEC std::vector<PCell<C,BT>*>
// store betti numbers for each frame
#define BETTI_STR std::map<BT, std::vector<BT>*>

// store (birth, death) intervals for each dimension
#define INTVEC std::vector<std::pair<BT,BT> >
#define INT_STR std::map<num, INTVEC*>

#define PCMAP std::map<Cell<C,BT>*, PCell<C,BT>*, CELLORD >

template <typename C = int, typename BT = int>
class PComplex
{
public:
	// complex data: persistent cells ordered by birth +  dimension
	PCVEC klist;
	INT_STR ints;
	BETTI_STR betti;

	// functions:
	PComplex()
	{
		ints.clear();
		betti.clear();
		klist.clear();
	}

	~PComplex()
	{
		destroyCells();
		destroyPersData();
	}


	friend std::ostream& operator << (std::ostream& out, const PComplex<C,BT>& toprint)
	{
		out<<"\n\nPCOMPLEX PRINTER...";
		for (num i = 0; i < (num) toprint.klist.size(); ++i)
		{
			out<<"\n   "<<*(toprint.klist.at(i))<<" bd "<<toprint.klist.at(i)->getBD();
		}
		return out;
	}

	num makeFromComplex(const Complex<C,BT>&, bool);
	bool makeBDChains (PCMAP&,bool);

	// Persistence Routines
	void showBetti(std::ostream&); // shows betti numbers for each frame to std out
	std::map<num, std::vector<std::pair<BT,BT> > > getInts() const; // returns persistence intervals
	void showInts(std::ostream&); // shows persistence intervals for all dimensions
	bool makeOutputFiles(const std::string&); // writes output to files
	void initPersData(const Complex<C,BT>&);
	void destroyCells();
	void destroyPersData();
	
	///PERSISTENCE ROUTINES
	/*num incrementBetti(const BT&, const BT&, const num&);
	void REMOVE_PIVOT_ROW(const PCell<C,BT>*, PCCHAIN&, bool);
	void COMPUTE_INTERVALS(const Complex<C,BT>&, bool, bool);*/
	num incrementBetti(const BT& from, const BT& to, const num& dim);
	void REMOVE_PIVOT_ROW(const PCell<C,BT>* cell, PCCHAIN& delta, bool trace);
	void COMPUTE_INTERVALS(const Complex<C,BT>& other, bool makegens, bool truncate);
};

#endif /* PCOMPLEX_H_ */
