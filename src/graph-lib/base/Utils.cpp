/// ***********************************************************
///	File	: Utils.cpp
///	About	: Useful methods to improve the execution of the application
///	Author	: mitxael@hotmail.it
/// License	: Copyright (C) 2018 Michael Vasquez Otazu, MIT License (MIT)
/// ******************************************************************

#include "Utils.h"				/// Main program

/// Definition of variables, to be used as global elsewhere through "extern"
	int verbosity;

void Utils::startExec()
{
	/// Read parameters
	int algorithm = 0;
	std::cout << "Select algorithm: \n1:Horton; 2:DePina; 3:Hybrid; 4:Amaldi; 5:Criccaldi; 6:Bron-Kerbosch (naive); 7:Bron-Kerbosch (tomita); 8:Bron-Kerbosch (eppstein); 9:PH (javaplex); 10:PH (dipha);";	///Dharwadker
	while (algorithm < 1 || algorithm > 11)			std::cin >> algorithm;
	
	std::cout << "Select verbosity (0:none; 1:normal; 2:detailed)";		//3:"exhaustive" is for debug only
	int verb = -1;
	while (verb < 0 || verb > 3)	std::cin >> verb;
	verbosity = verb;
	
	std::cout << "Select graph input (1:File; 2:Hypercube; 3:Random w/density 4: Euclidean)";
	int gInput = 0;
	while (gInput < 1 || gInput > 4)	std::cin >> gInput;

	boolean save2file = false;
	if(gInput > 1)
	{
		std::cout << "Save graph into file? (y/n)";
		char answer = '0';
		while (answer != 'y' && answer != 'n' && answer != '\n')	std::cin.get(answer) >> answer;
		if (answer == 'y') save2file = true;
	}

	/// Create Graph
	Graph G1(0);			/// Create an empty graph.
	switch (gInput) {
	case 1:					/// create a Graph from a selection menu (files in a directory)
	{
		G1.buildGraph();
		break;
	}
		/// create a Graph as a 2^n hypercube
		case 2:					
	{
		int n = 0;
		std::cout << "Enter number of vertices: (n>=3)";
		while (n < 3)	std::cin >> n;
		G1.buildGraph(save2file, n);
		int n_2 = (int)pow(2, n);
		break;
	}
		/// create a random Graph with given vertices and density (uniformly-distributed edges)
		case 3:					
	{
		std::cout << "Enter number of vertices: (n>=3)";
		int n = 0;
		while (n < 3)	std::cin >> n;
		double min = (double)(2.00000000 / (n - 1));					// minimum value required to create a connected graph
		std::cout << "Enter the density: (" + std::to_string(min) + " < density <= 1.0000)";
		double density = 0;
		while (density < min || density > 1)	std::cin >> density;
		G1.buildGraph(save2file, n, density*1.0000);
		break;
	}
	
		/// create an euclidean Graph with given vertices and edges (randomly added)
		case 4:					
	{
		std::cout << "Enter number of vertices: (n>=3)";
		int n = 0;
		while (n < 3)	std::cin >> n;
		int min = n;
		int max = (n)*(n - 1) / 2;
		std::cout << "Enter number of edges: ("<< std::to_string(min) << " < m < " << max << ")";
		int m = 0;
		while (m < 3 || m > max)	std::cin >> m;
		G1.buildGraph(save2file, n, m);
		break;
	}
	default:
	{
		std::cout << "Entered method is not valid!";
		break;
	}
	}

	/// Execute algorithm on Graph
	bool timer = true;		/// set the timer mode (true:high-resolution; false:standard)
	bool saveToFile = 0;	/// set the output mode (true:file; false:console)
	std::vector<std::pair<std::string, std::string>> algoParams = std::vector<std::pair<std::string, std::string>>();
	double Runtime = this->measurePerformance(G1, algorithm, verb, timer, saveToFile, algoParams);
	if (saveToFile != false)
		std::cout << std::endl << "**** CPU TIME: " << Runtime << "seconds ****" << std::endl << std::endl;
	
	//std::system("pause");
}

double Utils::measurePerformance(Graph &G1, int _algorithm, int verb, bool hiresClock, bool saveToFile, 
								std::vector<std::pair<std::string, std::string>> algoParams)
{
	/// if save2file is required, redirect console's output
	int old_stdout = (saveToFile == true)? outputToFile_Take() : 0;

	double RuntimeSec = 0;
	verbosity = verb;

	if (G1.V >= 3 && G1.E >= G1.V)
	{
		if ( (_algorithm > 4) || (_algorithm <= 4 && G1.isConnected("bfs")) )
		{
			std::cout << "***********************************************************************" << std::endl;
			std::cout << "The *" << G1.algorithmS[_algorithm-1] << "* algorithm is being executed for G: n = " << G1.V << " m = " << G1.E << " v = " << G1.E - G1.V + 1 << std::endl;
			std::cout << "(source = " << G1.source << ")" << std::endl;
			std::cout << "***********************************************************************" << std::endl;
			writeLog("***************** New Execution for <<" + G1.algorithmS[_algorithm-1] + ">> *********************");

			/// display graph
			if (verbosity > 0) G1.displayGraph();

			/// initialize Graph's ancillary elements
			G1.reset();

			/// activate execution measurement
			auto startTime = startChrono<std::chrono::high_resolution_clock::time_point>();

			/// Execute algorithm
			switch (_algorithm) {
			case 1:
			{
				Horton(G1);
				break;
			}
			case 2:
			{
				DePina(G1);
				break;
			}
			case 3:
			{
				Hybrid(G1);
				break;
			}
			case 4:
			{
				Amaldi(G1);
				break;
			}
			case 5: case 6: case 7:
			{
				/// Get parameters
				int variant = std::stoi(algoParams.at(0).second);				/// 0:Naive; 1:Tomita; 2:Eppstein
				int kMin = std::stoi(algoParams.at(1).second) + 1;				/// 1,2,3...

				/// Find cliques
				BronKerbosch(G1, variant, kMin);
				break;
			}
			case 8:
			{
				/// Get parameters
				int kMin = std::stoi(algoParams.at(1).second) + 1;				/// 1,2,3...

				/// Find cliques
				Criccaldi(G1);
				break;
			}
			//case 9:
			//{
				//Dharwadker(G1);
				//break;
			//}
			case 9: case 10: case 11: case 12: case 13:
			{
				/// Get parameters
				int phEngine = std::stoi(algoParams.at(0).second);				/// 0:Javaplex; 1:PHAT; 2:DiPHA; 3:Perseus; 4:Gudhi
				int complexType = std::stoi(algoParams.at(1).second);			/// 0:AllCliques; 1:Graph; 2:MaxCliques; 3:MinCycles
				bool onlyInfinite = (std::stoi(algoParams.at(2).second) == 0);	/// 0: Yes; 1:No
				int dimension = std::stoi(algoParams.at(3).second) - 1;			///-1:All; 0:dim-0; 1:dim-1
				bool saveCycles = (std::stoi(algoParams.at(4).second) == 0);	/// 0: Yes; 1:No

				/// Compute PH
				Homology::PersistentHomology(G1, phEngine, complexType, onlyInfinite, dimension, saveCycles);	
				break;
			}
			default:
			{
				std::cout << " Entered algorithm is not valid!";
				break;
			}
			}

			/// Deactivate execution measurement
			RuntimeSec = stopChrono(startTime);

			/// Display info for MCB algorithms
			if(_algorithm < 5)
			{
				/// verify the correctness of the MCB
				bool success = G1.verifyMCB();
				success? std::cout << "A consistent MCB was found!\n" : std::cout << "No consistent MCB was found.\n" << std::endl;

				/// Display MCB depending on verbosity
				if (success) double MCB_weight = G1.displayMCB(true);
			}

			/// Display algorithm timing
			std::cout << std::endl << "**** CPU TIME: " << RuntimeSec << " seconds ****" << std::endl << std::endl;
			//writeResults("**** CPU TIME: " + std::to_string(RuntimeSec) + " seconds ****", true);
		}
		else
		{
			std::cout << "Unable to calculate the MCB. The Graph is not connected." << std::endl;
			}
	}
	else
	{
		std::cout << "Unable to calculate the MCB. Likely causes: (1) Wrong format (2) Insufficient vertices (n<3) (3) Insufficient edges (m<n)" << std::endl;
	}

	/// If changed, redirect console's output back
	if (old_stdout > 0) outputToFile_Release(old_stdout);

	return RuntimeSec;
}

int Utils::outputToFile_Take()
{
	int _old_stdout = _dup(1);								/// _dup: save file descriptor (1==console) into variable
	std::string filename = setLogFilename();
	FILE *fp1;
	freopen_s(&fp1, filename.c_str(), "w", stdout);			/// open file and stream to it  (w=createFile; a=appendToFile)
	return _old_stdout;
}

int Utils::outputToFile_Release(int old_stdout)
{
	fclose(stdout);											/// fclose(fp1)????
	FILE *fp2 = _fdopen(old_stdout, "w");
	freopen_s(&fp2, "CONOUT$", "w", stdout);				/// open file and stream to it  (w=createFile; a=appendToFile)
	_dup2(old_stdout, 1);									/// _dup2: restore original file descriptor
	return 0;
}

std::string Utils::setLogFilename()
{
	int counter = 0;
	std::string directory = ".\\output\\";
	bool cd = std::experimental::filesystem::create_directory(directory);
	for (auto &f : std::experimental::filesystem::directory_iterator(directory))
	{
		if ( (f.path().filename().generic_string().find("result_", 0) == 0) && (f.path().extension() == ".txt") )
		{
			/// Get filename into a string-stream
			std::stringstream ss;
			ss << f;
			std::string filename = ss.str();

			/// Get the number from filename
			std::string::size_type filename_prefix = filename.find('_');
			std::string::size_type filename_suffix = filename.find('.', filename_prefix + 2);
			int n = std::stoi(filename.substr(filename_prefix + 1, filename_suffix - filename_prefix - 1));
			counter = (n > counter) ? n : counter;
		}
	}

	counter = counter + 1;
	std::stringstream new_filename;
	new_filename << ".\\output\\result_" << std::to_string(counter) << ".txt";
	return new_filename.str();
}

template<typename T>
void Utils::writeResults(const T &data, bool firstLine)
{
	/// Determine filename only the first time
	///std::string current_result = (firstLine == true)? getFilename() : a-global-variable;
	std::string current_result = ".\\output\\result_1.txt";

	std::ofstream fout(current_result);
	if (fout.is_open()) {
		fout << data;
		fout << "\n";
		fout.close();
	}
}

void Utils::terminateExecution()
{
	exit(0);
	//TerminateThread();
	//abort();
	//std::terminate();
}

std::wstring Utils::ConvertStringToWString(std::string s)
{
	///INLINE WAY #1: L"abcde..."
	///INLINE WAY #2: TEXT("abcde...")
	int len;
	int slength = (int)s.length() + 1;
	len = MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, 0, 0); 
	wchar_t* buf = new wchar_t[len];
	MultiByteToWideChar(CP_ACP, 0, s.c_str(), slength, buf, len);
	std::wstring res(buf);
	delete[] buf;
	return res;
}

std::string Utils::GetStringRegKey(std::string path, std::string value)
{
	///Convert STRING to LPCWTR
	std::wstring tmp1 = ConvertStringToWString(path);
	auto spath = tmp1.c_str();

	///Convert STRING to LPCWTR
	std::wstring tmp2 = ConvertStringToWString(value);
	auto svalue = tmp2.c_str();

	///Get REG_SZ
	///KEY_ALL_ACCESS | KEY_WOW64_64KEY
	HKEY hKey;
	TCHAR buffer[100];
	DWORD size = sizeof(buffer);
	long lError = RegOpenKeyEx(HKEY_LOCAL_MACHINE, spath, 0, KEY_READ, &hKey);
	long lStatus = RegQueryValueEx(hKey, svalue, 0, NULL, (BYTE*) buffer, &size);
	std::wstring wkey(buffer);
	std::string skey(wkey.begin(), wkey.end());
	RegCloseKey(hKey);

	if(lError == 0 && lStatus == 0)
		return skey;
	else
		return "";
}
