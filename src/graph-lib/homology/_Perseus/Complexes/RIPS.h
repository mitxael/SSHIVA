/*
 * RIPS.h
 *
 *  Created on: Dec 2, 2010
 *      Author: Vidit
 */

#ifndef RIPS_H_
#define RIPS_H_

#include <algorithm> // for set intersections!
#include <cmath> // for ceil and floor etc.

#include "SToplex.h"
// upper triangular matrix to store point-point neighbor relations
#define RROW std::vector<std::pair<num,BT> >
#define RIPMAT std::vector<RROW*>
// vector of vertices
#define VERTS std::vector<Point<PS>*>

// structure to store high dimensional (> 0) simplices
#define SLIST std::vector<Simplex<C,PS,BT>*>
#define INSRIP std::map<num, SLIST*>

// point-birth pairs from file:
# define PB_PAIR std::map<Point<PS>*, std::pair<double,BT>, PT_ORDER >


template <typename C=int, typename PS=double, typename BT=int>
class RIPS : public SToplex<C,PS,BT>
{
public:
	// DATA
	RIPMAT isnbr; // contains upper-triangular dist matrix info
	VERTS vert; // contains vertices

	std::vector<double> rads; // contains radius for each point
	std::vector<BT> birth; // vector of birth times in case of

	double epsilon, stepsize;
	num numsteps;

	// FUNCTIONS

	RIPS()
	{
		InitRIPS();
	}

	void InitRIPS()
	{
		epsilon = 1; stepsize = 1;
		numsteps = 0;
		isnbr.clear();
		vert.clear();
		rads.clear();
		// in case we begin to allocate points...
		srand(time(NULL));
	}

	~RIPS()
	{
		clear();
	}

	void allocateNbrMatrix();
	num addVertex(const BT&);
	bool addEdge(const num, const num, const BT&);

	BT getEdgeBirth(num, num) const;
	void clear(); // clean rips structure!
	std::pair<num,bool> makeFromFile(std::ifstream&); // make rips complex where births are inherited from vertices
	std::pair<num,bool> makeFromFile_GrowBalls(std::ifstream&,bool,bool,bool); // make rips complex by growing epsilon balls around points
	std::pair<num,bool> makeFromDistMatrixFile(std::ifstream&,bool,bool); // make rips complex from matrix of pairwise distances
	std::pair<num,bool> makeFromTimeSeriesFile(std::ifstream&); // make rips complex from matrix of pairwise distances

	num witness(double); // prune vertex via witnesses. call before building nbrmatrix!
	void getLowerNbrs(num, std::vector<num>&) const; // get lower order neighbors!
	void getLowerNbrs_GrowBalls(num, std::vector<num>&, std::vector<BT>&, bool) const;
	void incRIPS(const num); // builds RIPS complex incrementally
	void incRIPS_GrowBalls(const num, const BT&);
	void makeRandom_EdgeShuffle(const num, bool);

	bool AddRIPSCofaces(Simplex<C,PS,BT>*, const std::vector<num>&, num, INSRIP&, num);
	bool AddRIPSCofaces_GrowBalls(Simplex<C,PS,BT>*, const std::vector<num>&, const std::vector<BT>&, num, INSRIP&);

	void makeNbrMatrix(); // populates upper triangular matrix of point neighbor relations
	void makeNbrMatrix_GrowBalls();

	void clearPoints(); // removes memory allocated to vertices

	void showPoints(const std::vector<num>&) const; // debug help for nbr lists!
	void storeCoface(const Simplex<C,PS,BT>*, INSRIP&) const;
	void insertToToplex(INSRIP&);
	void ComputePersistence(num,std::string,bool);
	void ComputePersistence_GrowBalls(num,std::string);

	void ComputePersistence(num ptdim, std::map<num, std::vector<std::pair<BT,BT> > >&);
	void showNbrMatrix(std::ostream& out) const;
};

#endif /* RIPS_H_ */
