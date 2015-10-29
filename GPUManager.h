#ifndef GPU_MANAGER_H
#define GPU_MANAGER_H

#include <cuda_runtime_api.h>
#include <cuda.h>
#include <cmath>
#include "cross_kernal.cuh"
#include "Input.h"



class GpuManager{
	private:
		size_t available_memory;
		size_t total_memory;
		//GPU pointers
		double* g_freq;
		double* g_intens;
		double* g_energies;
		double* g_nu;
		double * g_aif;
		double * g_gamma;
		double * g_n;
		int* g_gns;
		int N_trans;
		bool alloc;
		ProfileType profile;
		void ExecuteVoigtCrossSection(int N, int N_ener,int start_idx);
		void ExecuteVoigtCrossSectionBlock(int N, int N_ener,int start_idx);
		void ExecuteGaussianCrossSection(int N, int N_ener,int start_idx);	
		void ExecuteDopplerCrossSection(int N, int N_ener,int start_idx);	

	public:
		GpuManager(ProfileType pProfile);
		~GpuManager();
		void InitializeVectors(int Npoints);
		void InitializeConstants(double half_width,double temperature, double partition,double dfreq,double meanmass,double pressure=10.0,double ref_temp=296.0);
		void TransferVectors(size_t Nener,double* h_energies, double* h_nu, double* h_aif,int* h_gns,double* h_gamma=NULL,double * h_n=NULL);
		void TransferFreq(double* h_freq,double* h_intens,int N);
		void ExecuteCrossSection(int N, int N_ener,int start_idx);
		void TransferResults(double* h_freq,double* h_intens,int N);
		void Cleanup();
		void TrackMemory(size_t bytes);
		void FreeMemory(size_t bytes);
		size_t GetAvailableMemory();
		int GetNtrans(){ return N_trans;}
		

};

#endif
