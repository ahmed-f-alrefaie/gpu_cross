#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include "Util.h"
#pragma once
using namespace std;

enum ProfileType{
	GAUSSIAN,
	DOPPLER,
	VOIGT,
	PSUEDO_VOIGT
	
};

enum FileType{
	HITRAN_TYPE,
	EXOMOL_TYPE
};

class Input{
	private:
		double half_width;
		double mean_mass;
		double lorentz;
		double temperature;
		double partition;
		double pressure;
		double nu_start;
		double nu_end;
		int Npoints;
		int num_files;
		double gamma_air;
		double n_air;
		double hitran_mixture_air;
		int max_points;
		FileType which_file;
		vector<string> trans_files;
		vector<string> broadener_files;
		string state_file;
		string hitran_file;
		ProfileType StringToProfile(string & pString);
		//istream file_stream;
		
		
	public:
		Input();
		~Input();
		ProfileType profile;
		void ReadStream();
		void ReadInput();
		void ReadInputII();
		double GetHalfWidth(){return half_width;};
		double GetMeanMass(){return mean_mass;};
		double GetTemperature(){return temperature;};
		double GetNuStart(){return nu_start;};
		double GetNuEnd(){return nu_end;};
		double GetPartition(){return partition;};
		double GetPressure(){return pressure;}
		double GetHitranMixture(){return hitran_mixture_air;}
		double GetMaxPoints(){return max_points;}
		FileType GetFileType(){return which_file;}

		int GetNpoints(){return Npoints;};
		string GetStateFile(){return state_file;};
		vector<string> GetTransFiles(){return trans_files;};
		vector<string> GetBroadeners(){return broadener_files;}
		








};
