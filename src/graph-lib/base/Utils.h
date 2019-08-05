/// ***********************************************************
///	File	: Utils.h
///	About	: Header of Utils class
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
/// ******************************************************************


#ifndef UTILS_H
#define UTILS_H
#pragma once

#include "Graph.h"					/// Graph's header file
#include "../cycles/Horton.h"					/// Horton's header file
#include "../cycles/DePina.h"					/// DePina's header file
#include "../cycles/Hybrid.h"
#include "../cycles/Amaldi.h"					/// Amaldi's header file
#include "../cliques/BronKerbosch.h"			/// BronKerbosch's header file
#include "../cliques/Criccaldi.h"				/// Criccaldi's header file
#include "../cliques/Dharwadker.h"				/// Dharwadker's header file
#include "../homology/Homology.h"				/// Homology header file

#include <iostream>					/// std::cin std::cout
#include <chrono>					/// std::chrono
#include <ctime>					/// std::clock_t
#include <io.h>						/// input output
#include <filesystem>				/// iterate into the input directoty
//#include <thread>					/// Execute algorithms in a new thread

#define LOG(msg) if( verbose ) std::cout << msg << std::endl;

class Utils
{
public:

	/// Definition of methods
	void startExec();																/// Prepare data and execute algorithm
	double measurePerformance(Graph &G1, int algorithm, int verbosity, bool hiresClock, 
		bool saveToFile, std::vector<std::pair<std::string, std::string>> algoParams); /// Measure the execution performance of an algorithm
	std::string setLogFilename();													/// Determines the name of the log file
	int outputToFile_Take();														/// Take console, and redirect its output to a file
	int outputToFile_Release(int old_stdout);										/// Release console, and redirect its output back
	template<typename T> static void writeResults(const T &data, bool firstLine);	/// Write results on the results.txt file
	void terminateExecution();														/// Force the execution to terminate

	/// Definition of static methods
	template<typename T> static void writeLog(const T &data)
	{
		///out:open for writing	| app:append at end
		std::ofstream flog;
		flog.open(".\\graph-suite.log", std::ofstream::out | std::ofstream::app);
		char datetime[100]; time_t now = time(0); struct tm buf; localtime_s(&buf, &now); strftime(datetime, 100, "%Y-%m-%d %H:%M:%S.000", &buf);
		flog << datetime << "\t" << data << std::endl;
		flog.close();
	}

	template<typename T> static T startChrono();
	template<>
	static std::chrono::high_resolution_clock::time_point startChrono()
	{
		std::chrono::high_resolution_clock::time_point startTime;
		startTime = std::chrono::high_resolution_clock::now();
		return startTime;
	}
	template<>
	static std::clock_t startChrono()
	{
		std::clock_t startTime;
		startTime = std::clock();
		return startTime;
	}

	template<typename T> static double stopChrono(T startTime);
	template<>
	static double stopChrono(std::chrono::high_resolution_clock::time_point startTime)
	{
		std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
		/// calculate in microseconds
		auto duration_chrono_us = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
		/// convert into seconds
		double RuntimeSec = duration_chrono_us / double(1000000.00000000);
		return RuntimeSec;
	}
	template<>
	static double stopChrono(std::clock_t startTime)
	{
		std::clock_t endTime = std::clock();
		/// calculate in microseconds:
		auto duration_clock_ms = (endTime - startTime) / CLOCKS_PER_SEC;
		/// convert into seconds
		double RuntimeSec = duration_clock_ms / double(1000000.00000000);
		return RuntimeSec;
	}
	
	static std::wstring ConvertStringToWString(std::string s);
	
	static std::string GetStringRegKey(std::string path, std::string value);

	/// Definition of friends
	// friend ... ...(... ...);
};
#endif