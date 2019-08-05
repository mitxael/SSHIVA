/// ******************************************************************
///	File	: Dharwadker.cpp
///	About	: Dharwadker's method for finding the Maximal Cliques
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
///
/// Ref.	:
/// -  New polynomial-time algorithm for finding maximal cliques in graphs (A. Dharwadker) (2006)
/// - (http://www.dharwadker.org/clique/clique.pdf)
/// ******************************************************************

#include "Dharwadker.h"				/// Cliques's header file

bool removable(std::vector<int> neighbor, std::vector<int> cover);
int max_removable(std::vector<std::vector<int> > neighbors, std::vector<int> cover);
std::vector<int> procedure_1(std::vector<std::vector<int> > neighbors, std::vector<int> cover);
std::vector<int> procedure_2(std::vector<std::vector<int> > neighbors, std::vector<int> cover, int k);
int cover_size(std::vector<int> cover);
typedef std::vector<std::vector<int>> _graph;

std::tuple<int, int, int, _graph, _graph> readGraph(std::string filePath)
{
	/// Read Graph (note we work with the complement of the input Graph) 
	std::ifstream infile(filePath);
	int n, i, j, k, K, edge = 0;
	infile >> n;
	std::vector<std::vector<int>> graph;

	/// Read adjacent list
	for (i = 0; i<n; i++)
	{
		std::vector<int> row;
		for (j = 0; j<n; j++)
			row.push_back(1);
			graph.push_back(row);
		}
		int m, r, c;
		double w;
		infile >> m;
		for (int z = 0; z < m; z++) {
			infile >> r >> c >> w;
			graph[c][r] = 0;
			graph[r][c] = 0;
	}		

	/// Find Neighbors 
	std::vector<std::vector<int>> neighbors;
	for (i = 0; i<graph.size(); i++)
	{
		std::vector<int> neighbor;
		for (j = 0; j<graph[i].size(); j++)
			if (graph[i][j] == 1)
				neighbor.push_back(j);
		neighbors.push_back(neighbor);
	}

	/// Set maximum size of Clique wanted 
	K = 3;
	k = n - K;

	return std::make_tuple(n , k , K, graph, neighbors);
}

int Dharwadker(Graph& G)
{
	std::string inputFile = G.source;

	/// Check input
	if(!inputFile.empty())
	{
		inputFile = ".\\input\\" + inputFile;
	}
	else
	{
		std::cout << "Wrong input file!" << std::endl;
		return -1;
	}

	/// Create graph
	int i, j, p, q, r, s, min, counter = 0;
	std::tuple<int, int, int, _graph, _graph> newGraph = readGraph(inputFile);
	int n = std::get<0>(newGraph);
	int k = std::get<1>(newGraph);
	int K = std::get<2>(newGraph);
	std::vector<std::vector<int>> graph = std::get<3>(newGraph);
	std::vector<std::vector<int>> neighbors = std::get<4>(newGraph);
	
	/// Find Cliques 
	bool found = false;

	std::cout << "Finding Cliques..." << std::endl;
	min = n + 1;
	std::vector<std::vector<int> > covers;
	std::vector<int> allcover;

	for (i = 0; i < graph.size(); i++)
		allcover.push_back(1);
	
	for (i = 0; i < allcover.size(); i++)
	{
		if (found)
			break;
		
		counter++;
		
		std::cout << counter << ". ";
		
		std::vector<int> cover = allcover;
		cover[i] = 0;
		cover = procedure_1(neighbors, cover);
		s = cover_size(cover);
		
		if (s<min)
			min = s;

		if (s <= k)
		{
			std::cout << "Clique (" << n - s << "): ";
			for (j = 0; j<cover.size(); j++)
				if (cover[j] == 0)
					std::cout << j + 1 << " ";
			std::cout << std::endl;
			std::cout << "Clique Size: " << n - s << std::endl;
			covers.push_back(cover);
			found = true;
			break;
		}
		
		for (j = 0; j<n - k; j++)
			cover = procedure_2(neighbors, cover, j);
		s = cover_size(cover);
		
		if (s<min)
			min = s;
		
		std::cout << "Clique (" << n - s << "): ";
		
		for (j = 0; j<cover.size(); j++)
			if (cover[j] == 0)
				std::cout << j + 1 << " ";
		std::cout << std::endl;

		std::cout << "Clique Size: " << n - s << std::endl;
		
		covers.push_back(cover);
		
		if (s <= k)
		{
			found = true;
			break;
		}
	}

	/// Pairwise Intersections 
	for (p = 0; p < covers.size(); p++)
	{
		if (found) break;
		for (q = p + 1; q < covers.size(); q++)
		{
			if (found)
				break;
			
			counter++;
			
			std::cout << counter << ". ";
			std::vector<int> cover = allcover;
			
			for (r = 0; r<cover.size(); r++)
				if (covers[p][r] == 0 && covers[q][r] == 0) cover[r] = 0;
			cover = procedure_1(neighbors, cover);
			s = cover_size(cover);
			
			if (s < min)
				min = s;

			if (s <= k)
			{
				std::cout << "Clique (" << n - s << "): ";
				for (j = 0; j<cover.size(); j++) 
					if (cover[j] == 0) 
						std::cout << j + 1 << " ";
				std::cout << std::endl;
				std::cout << "Clique Size: " << n - s << std::endl;
				found = true;
				break;
			}

			for (j = 0; j<k; j++)
				cover = procedure_2(neighbors, cover, j);
			s = cover_size(cover);

			if (s < min)
				min = s;
			
			std::cout << "Clique (" << n - s << "): ";

			for (j = 0; j<cover.size(); j++)
				if (cover[j] == 0)
					std::cout << j + 1 << " ";
			std::cout << std::endl;
			
			std::cout << "Clique Size: " << n - s << std::endl;

			if (s <= k)
			{
				found = true;
				break;
			}
		}
	}

	std::cout << "Maximum Clique size found is " << n - min << "." << std::endl;

	system("PAUSE");
	
	return 0;
}

bool removable(std::vector<int> neighbor, std::vector<int> cover)
{
	bool check = true;
	for (int i = 0; i<neighbor.size(); i++)
		if (cover[neighbor[i]] == 0)
		{
			check = false;
			break;
		}
	return check;
}

int max_removable(std::vector<std::vector<int> > neighbors, std::vector<int> cover)
{
	int r = -1, max = -1;
	for (int i = 0; i<cover.size(); i++)
	{
		if (cover[i] == 1 && removable(neighbors[i], cover) == true)
		{
			std::vector<int> temp_cover = cover;
			temp_cover[i] = 0;
			int sum = 0;
			for (int j = 0; j<temp_cover.size(); j++)
				if (temp_cover[j] == 1 && removable(neighbors[j], temp_cover) == true)
					sum++;
			if (sum>max)
			{
				max = sum;
				r = i;
			}
		}
	}
	return r;
}

std::vector<int> procedure_1(std::vector<std::vector<int> > neighbors, std::vector<int> cover)
{
	std::vector<int> temp_cover = cover;
	int r = 0;
	while (r != -1)
	{
		r = max_removable(neighbors, temp_cover);
		if (r != -1) temp_cover[r] = 0;
	}
	return temp_cover;
}

std::vector<int> procedure_2(std::vector<std::vector<int> > neighbors, std::vector<int> cover, int k)
{
	int count = 0;
	std::vector<int> temp_cover = cover;
	int i;
	for (i = 0; i<temp_cover.size(); i++)
	{
		if (temp_cover[i] == 1)
		{
			int sum = 0, index = 0;
			for (int j = 0; j<neighbors[i].size(); j++)
				if (temp_cover[neighbors[i][j]] == 0) { index = j; sum++; }
			if (sum == 1 && cover[neighbors[i][index]] == 0)
			{
				temp_cover[neighbors[i][index]] = 1;
				temp_cover[i] = 0;
				temp_cover = procedure_1(neighbors, temp_cover);
				count++;
			}
			if (count>k) break;
		}
	}
	return temp_cover;
}

int cover_size(std::vector<int> cover)
{
	int count = 0;
	for (int i = 0; i<cover.size(); i++)
		if (cover[i] == 1) count++;
	return count;
}