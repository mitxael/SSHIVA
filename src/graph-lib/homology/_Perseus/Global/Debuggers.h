/*
 * Debuggers.h
 *
 *  Created on: Dec 31, 2010
 *      Author: Vidit
 */

#ifndef DEBUGGERS_H_
#define DEBUGGERS_H_

#include<ostream>
#include<deque>
#include<vector>
#include <iostream>

//// prints vector of stuff
//template <typename data>
//void vprint(const vector<data>& toprint, ostream& out = cout)
//{
//	typename vector<data>::const_iterator p;
//	out<<"< ";
//	for (p = toprint.begin(); p != toprint.end(); ++p)
//	{
//		out<<*p<<" ";
//	}
//	out<<">";
//}

// prints container of stuff
template <typename container>
void c_print (const container& toprint, std::ostream& out = std::cout)
{
	typename container::const_iterator cit;
	//out<<"< ";
	for (cit = toprint.begin(); cit != toprint.end(); ++cit)
	{
		out<<*cit<<" ";
	}
	//out<<">";
}

// prints container of pointers to stuff
template <typename container>
void cptr_print (const container& toprint, std::ostream& out = std::cout)
{
	typename container::const_iterator cit;
	out<<"< ";
	for (cit = toprint.begin(); cit != toprint.end(); ++cit)
	{
		out<<*(*cit)<<" ";
	}
	out<<">";
}


// displays queue for debugging purposes!
template <typename data_ptr>
void printQ(std::deque<data_ptr>& myq, std::ostream& out = std::cout)
{
	out<<"\n*******Q*******\n";
	typename std::deque<data_ptr>::const_iterator qiter;
	for (qiter = myq.begin(); qiter != myq.end(); ++qiter)
	{
		out<<*qiter<<" ";
	}
	out<<'\n';
}

#endif /* DEBUGGERS_H_ */
