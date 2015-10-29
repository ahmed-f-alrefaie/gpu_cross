#include "MultiGPUManager.h"
#include "StateReader.h"
#include "Input.h"
#include "Timer.h"

#pragma once


class BaseProfile{
	protected:
		MultiGpuManager* gpu_manager;
		StateReader* state_reader;
		Input* input;
		ProfileType profile;
		double* freq;
		double* intens;
		double dfreq;
		
		double partition;

		int ib,ie;
		double* h_energies; 
		double* h_nu;
		double * h_aif;
		int* h_gns;
		double* h_gammaL;
		double* h_n;

		int Npoints;
		int num_trans_fit;
		int Ntrans;
		int points;
		int old_Ntrans;
		double start_nu,end_nu,min_nu,max_nu,max_hw;
		vector<string> transition_files;
		BaseProfile(Input* pInput);
		~BaseProfile();
		virtual void PerformCrossSection(double HW);
		
		virtual void InitializeProfile()=0;
	public:
		void Initialize();
		virtual void ExecuteCrossSection()=0;
		void OutputProfile();
};
