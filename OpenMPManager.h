#include <cmath>
#include <thread>
#include "BaseManager.h"
#include <vector>
#include "defines.h"
#pragma once

//This is a lie it is not openmp but uses threads
class OpenMPManager : public BaseManager{

private:
		bool* threads_done;
		int total_threads;
		std::vector<std::thread> workers;
		double** t_intens; 
		double m_half_width;
		double m_beta;
		double m_partition;
		double m_ref_temp;
		double m_dfreq;
		double m_temperature;
		double m_mean_mass;
		double m_pressure;	
		double m_cmcoef;
		double* g_freq;
		double* g_intens;
		double* g_energies;
		double* g_nu;
		double * g_aif;
		double * g_gamma;
		double * g_n;
		int* g_gns;	
		double start_nu;	
		void JoinAllThreads();
		void ComputeVoigt(int id,int Npoints,int Nener,int start_idx);		

public:
		OpenMPManager(ProfileType pprofile, int num_threads,size_t total_memory);
		void InitializeVectors(int Npoints);
		void InitializeConstants(double half_width,double temperature, double partition,double dfreq,double meanmass,double pressure=10.0,double ref_temp=296.0);
		void TransferVectors(size_t Nener,double* h_energies, double* h_nu, double* h_aif,int* h_gns,double* h_gamma=NULL,double * h_n=NULL);
		void TransferFreq(double* h_freq,double* h_intens,int N);
		void ExecuteCrossSection(int N, int N_ener,int start_idx);
		void TransferResults(double* h_freq,double* h_intens,int N);
		void Cleanup();
		bool ReadyForWork();

		

};
