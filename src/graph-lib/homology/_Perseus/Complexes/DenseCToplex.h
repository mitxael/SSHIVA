/*
 * DenseCToplex.h
 *
 *  Attempt at an efficient way to compute persistence intervals on an integer cubical grid
 */

#ifndef DENSECTOPLEX_H_
#define DENSECTOPLEX_H_

#include "../Cells/Cell.h"
#include "../Global/Global.h"
#include "../DebugV.h"
#include "Morse.h"

// the cubes are stored as a 2d array. the second index is the lexico order # of the
// (lowest, leftest...) "anchor" point of the cube, and the first one is the lexico
// ordering number of the "addin" vector
# define CVEC std::vector<Cell<C,BT>*>
# define CGRID std::vector<CVEC*>

// C is the coefficient ring over which we build chains, PS is the underlying point space
// to describe the vertices of the cubes, eg int if all vertices are integers etc. BT is
// the space from which persistence parameters are chosen, the birth times.

template <typename C = int, typename BT = int>
class DenseCToplex
{
public:
	CGRID cube; // 2d vector representing all cubes, see #define above
	std::vector<num> extent; // vector of size topdim, extent[i] is "grid size along dim i"
	bool faced; // have faces been made yet?

// functions:

	// Constructor
	DenseCToplex()
	{
		cube.clear();
		extent.clear();
		faced = false;
		//Init(ext);
	}

	// Destructor
	~DenseCToplex()
	{
		Destroy();
	}

	// helpful printer!
	friend std::ostream& operator << (std::ostream& out, const DenseCToplex<C,BT>& toprint)
	{
		num addnum = toprint.cube.size();

		std::cout<<"\n\n*********DENSE CTOP PRINTER***********\n";
		std::vector<num> advec, anch;
		CVEC* curvec;

		for (num i=0; i<addnum; i++)
		{
			// get and print Addin vector!
			if(!toprint.getAddin(i,advec)) continue;
			c_print<std::vector<num> >(advec, out); out<<"\n\n";
			// now get the size of the cube list with this addin vector
			curvec = toprint.cube.at(i);
			if (curvec == NULL) continue;
			for (num j=0; j<(num)curvec->size();j++)
			{
				if (curvec->at(j) == NULL) continue;
				// extract anchor point!
				if(toprint.getAnchor(j,i,anch))
				{
					out<<"     ";
					c_print<std::vector<num> >(anch); out<<" birth "<<curvec->at(j)->birth;
					out<<"\n";
					//cin.get();
				}
			}
		}
		return out;
	}

	// allocator and deallocator
	void Init(const std::vector<num>&); // initialize to given extent vector
	void Destroy(); // deallocates memory

	// makes faces
	bool makeFaces(bool);
	void makeFaceLinks(Cell<C,BT>*, const std::vector<num>&, const std::vector<num>&);

	void addTopCube(const std::vector<num>&, const BT&);

	// addin & anchor getters and setters
	bool getAddin(num, std::vector<num>&) const; // recover addin vector from lexico position
	bool getAnchor(const num, const num, std::vector<num>&) const; // recover anchor point's coordinates from lexico pos + addin vec
	bool getAnchor(const num, const std::vector<num>&, std::vector<num>&) const; // recover anchor point's coordinates from lexico pos + addin vec
	num getAddinPos(const std::vector<num>&, std::vector<num>&) const; // lexico position of input addin vector
	num getAddinPosFast(const std::vector<num>& advec) const;
	num getAnchPos(const std::vector<num>&, const num) const; // lexico position of anchor point.
	void getExtForAddin(const std::vector<num>&, std::vector<num>&) const; // get extent vector for given addin vector
	num getAnchPos(const std::vector<num>&, const std::vector<num>&) const;
	bool getAnchorFast(num,  const std::vector<num>&, std::vector<num>&) const; // recover anchor point from lexico position
	num getAnchPosFast(const std::vector<num>&, const std::vector<num>&) const;
	// file input, complex output
	std::pair<num,bool> makeFromFile(std::ifstream&);
	void writeComplex(MComplex<C,BT>&, bool);
	bool quickWriteComplex(MComplex<C,BT>&, bool);

	void ComputePersistence(std::string);
};

#endif /* DENSECTOPLEX_H_ */
