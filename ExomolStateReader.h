#include "StateReader.h"
#include "exomol_objects.h"
#include "exomol_functions.h"
#include "Util.h"
#include <fstream>
#include <string>
#include <vector>
#include <map>
//#include <bzlib.h>
#pragma once



struct broadData{
	double gamma;
	double n;
};

struct BroadInfo{
	std::map<std::string,broadData> broad_map;
	std::vector<std::string> code_list;
	double fraction;
	double default_gamma;
	double default_n;
};
class ExomolStateReader : public StateReader{
	private:
		FILE* trans_file;
		//std::ifstream stream;
		exomol_states exomol;
		BroadInfo broad_info;
		//bz2_stream zip_stream;
		bool use_broadeners;
		bool use_compressed;
		double m_pressure;
		bool GetGammaN(state* state_i,state* state_f,double & gamma, double & n);
		void InitializeBroadener(std::string filename);
		void Decompress();
		void OpenBz2(std::string pFilename);
		void CloseBz2();
	public:
		ExomolStateReader(std::string pFilename="",double partition=-1.0,double pressure=1.0,std::string broadeners="");
		~ExomolStateReader();
		bool OpenFile(std::string pFilename);
		bool CloseFile();
		bool ReadNextState(double & nu,double & gns,double & e_i, double & aif, double & gam,double & n);
		bool ReadNextState(double & nu,double & gns,double & e_i, double & aif);
		double ComputePartition(double temperature);
};
