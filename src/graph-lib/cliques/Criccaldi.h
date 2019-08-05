/// ******************************************************************
///	File	: Criccaldi.h
///	About	: Header of Criccaldi class
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// - Efficient Deterministic Algorithms for Finding a Minimal Cycle Basis (E. Amaldi, C. Iuliano, R. Rizzi) (2010)
/// - (https://link.springer.com/chapter/10.1007/978-3-642-13036-6_30)
/// ******************************************************************

#ifndef CRICCALDI_H
#define CRICCALDI_H
#pragma once

#include "../cycles/Amaldi.h"
//#include "BronKerbosch.h"

	/// Declaration of methods
	void Criccaldi(Graph& G1);											/// Cliques from Amaldi
	void maximalCliques(Graph &G1);										/// Find maximal cliques among cliques

	void allCliquesBF(Graph &G1, int kmax);

#endif