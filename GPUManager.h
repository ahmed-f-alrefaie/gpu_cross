#ifndef GPU_MANAGER_H
#define GPU_MANAGER_H
#include "BaseManager.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <cmath>
#include "cross_kernal.cuh"
#include "Input.h"



class GpuManager : public BaseManager{
	private:
		int gpu_id;
		double* g_freq;
		double* g_intens;
		double* g_energies;
		double* g_nu;
		double * g_aif;
		double * g_gamma;
		double * g_n;
		int* g_gns;
		//GPU pointers
		bool alloc;
		void ExecuteVoigtCrossSection(int N, int N_ener,int start_idx);
		void ExecuteVoigtCrossSectionBlock(int N, int N_ener,int start_idx);
		void ExecuteGaussianCrossSection(int N, int N_ener,int start_idx);	
		void ExecuteDopplerCrossSection(int N, int N_ener,int start_idx);	

	public:
		GpuManager(ProfileType pProfile,int gpu_id=0);
		~GpuManager();
		void InitializeVectors(int Npoints);
		void InitializeConstants(double half_width,double temperature, double partition,double dfreq,double meanmass,double pressure=10.0,double ref_temp=296.0);
		void TransferVectors(size_t Nener,double* h_energies, double* h_nu, double* h_aif,int* h_gns,double* h_gamma=NULL,double * h_n=NULL);
		void TransferFreq(double* h_freq,double* h_intens,int N);
		void ExecuteCrossSection(int N, int N_ener,int start_idx);
		void TransferResults(double* h_freq,double* h_intens,int N);
		void Cleanup();
		bool ReadyForWork();
		

};

#endif
