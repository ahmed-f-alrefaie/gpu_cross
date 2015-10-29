#include "GPUManager.h"
#include "Timer.h"
#include "cuda_utils.cuh"

GpuManager::GpuManager(ProfileType pProfile,int pgpu_id) : alloc(false),g_freq(NULL),g_intens(NULL),g_energies(NULL),g_nu(NULL),g_aif(NULL),g_gns(NULL),g_gamma(NULL),g_n(NULL), gpu_id(pgpu_id){
	cudaSetDevice(gpu_id);
	cudaFree(0);
	//Get device properties
	cudaDeviceProp devProp;
        cudaGetDeviceProperties(&devProp, 0);
	//Now we compute total bytes
	available_memory = devProp.totalGlobalMem;
	total_memory = devProp.totalGlobalMem;
	available_memory -= 300000000l;
	//available_memory /=6l;	
	CheckCudaError("Get Device Properties");
	profile = pProfile;

}
GpuManager::~GpuManager(){
	Cleanup();
	
}


void GpuManager::InitializeConstants(double half_width,double temperature, double partition,double dfreq,double mean_mass,double pressure,double ref_temp){
	cudaSetDevice(gpu_id);
	cross_section_data cross_constants;

	double planck =6.6260693e-27;
	double vellgt=2.99792458e10;
	double boltz=1.380658e-16;
	double pi=std::acos(-1.0);
	printf ("HW = %12.6f temp=%12.6f partition=%12.6f\n",half_width,temperature, partition);
	printf("Doppler constant = %14.4E\n",sqrt(2.0*boltz*temperature/mean_mass)/vellgt);
	cross_constants.halfwidth = half_width;
	cross_constants.beta = planck*vellgt/(boltz*temperature);
	cross_constants.partition = partition;
	cross_constants.ref_temp= ref_temp;
	//cross_constants.ln2pi= sqrt(cross_constants.ln2/pi);
	cross_constants.cmcoef=1.0/(8.0*pi*vellgt);
	cross_constants.dfreq = dfreq;
	cross_constants.temperature = temperature;
	cross_constants.mean_mass = mean_mass;
	cross_constants.pressure = pressure;
	copy_intensity_info(&cross_constants);
	CheckCudaError("Copy Constants");

}

void GpuManager::InitializeVectors(int Npoints){
	cudaSetDevice(gpu_id);
	cudaMalloc((void**)&g_freq,sizeof(double)*size_t(Npoints));
	TrackMemory(sizeof(double)*size_t(Npoints));
	cudaMalloc((void**)&g_intens,sizeof(double)*size_t(Npoints));
	TrackMemory(sizeof(double)*size_t(Npoints));
	if(profile==VOIGT)
		N_trans = available_memory/(sizeof(double)+sizeof(double)+sizeof(double)+sizeof(int)+sizeof(double) + sizeof(double));
	else
		N_trans = available_memory/(sizeof(double)+sizeof(double)+sizeof(double)+sizeof(int));
	printf("Number of transitions in the GPU = %d\n",N_trans);
	cudaMalloc((void**)&g_energies,sizeof(double)*size_t(N_trans));
	TrackMemory(sizeof(double)*size_t(N_trans));
	cudaMalloc((void**)&g_nu,sizeof(double)*size_t(N_trans));
	TrackMemory(sizeof(double)*size_t(N_trans));
	cudaMalloc((void**)&g_aif,sizeof(double)*size_t(N_trans));
	TrackMemory(sizeof(double)*size_t(N_trans));
	cudaMalloc((void**)&g_gns,sizeof(int)*size_t(N_trans));
	TrackMemory(sizeof(int)*size_t(N_trans));
	if(profile==VOIGT){
		printf("Allocating Voigt stuff");
		cudaMalloc((void**)&g_gamma,sizeof(double)*size_t(N_trans));
		TrackMemory(sizeof(double)*size_t(N_trans));
		cudaMalloc((void**)&g_n,sizeof(double)*size_t(N_trans));
		TrackMemory(sizeof(double)*size_t(N_trans));
	}
	CheckCudaError("Initialize Vectors");
}

void GpuManager::TransferFreq(double* h_freq,double* h_intens,int N){
	cudaSetDevice(gpu_id);
	cudaMemcpy(g_freq,h_freq,sizeof(double)*size_t(N),cudaMemcpyHostToDevice);
	cudaMemcpy(g_intens,h_intens,sizeof(double)*size_t(N),cudaMemcpyHostToDevice);


	cudaDeviceSynchronize();
	CheckCudaError("Copy");
}
void GpuManager::TransferVectors(size_t Nener,double* h_energies, double* h_nu, double* h_aif,int* h_gns,double* h_gamma,double * h_n){
	cudaSetDevice(gpu_id);
	cudaDeviceSynchronize(); // Synchronize at the beginning
	CheckCudaError("Kernal");
	Timer::getInstance().StartTimer("MemCpy");
	printf("Performing Copy, Ntransitions = %d\n",Nener);
	cudaMemcpy(g_energies,h_energies,sizeof(double)*size_t(Nener),cudaMemcpyHostToDevice);
	cudaMemcpy(g_nu,h_nu,sizeof(double)*size_t(Nener),cudaMemcpyHostToDevice);
	cudaMemcpy(g_aif,h_aif,sizeof(double)*size_t(Nener),cudaMemcpyHostToDevice);
	cudaMemcpy(g_gns,h_gns,sizeof(int)*size_t(Nener),cudaMemcpyHostToDevice);
	if(h_gamma != NULL && h_n != NULL){
		printf("Transferred gammas\n");
		cudaMemcpy(g_gamma,h_gamma,sizeof(double)*size_t(Nener),cudaMemcpyHostToDevice);
		cudaMemcpy(g_n,h_n,sizeof(double)*size_t(Nener),cudaMemcpyHostToDevice);
	}
	cudaDeviceSynchronize();
	CheckCudaError("Copy");
	Timer::getInstance().EndTimer("MemCpy");
}

void GpuManager::ExecuteCrossSection(int N, int N_ener,int start_idx){
	cudaSetDevice(gpu_id);
	if(profile == GAUSSIAN)
		execute_two_step_kernal_block(g_freq, g_intens,g_energies,g_nu,g_gns,g_aif,N,N_ener,start_idx);
	else if(profile == DOPPLER)
		//printf("Not implemented\n");
		ExecuteDopplerCrossSection(N, N_ener,start_idx);
	else if(profile == VOIGT){
		//if(N > 100000)
		//	ExecuteVoigtCrossSection(N, N_ener,start_idx);
		//else
			ExecuteVoigtCrossSectionBlock(N, N_ener,start_idx);
	}
}

void GpuManager::ExecuteGaussianCrossSection(int N, int N_ener,int start_idx){
		cudaSetDevice(gpu_id);		
		execute_two_step_kernal_block(g_freq, g_intens,g_energies,g_nu,g_gns,g_aif,N,N_ener,start_idx);
					//device_compute_cross_section_abscoef<<<gridSize,blockSize>>>(g_freq, g_cs,g_energies,g_gns,g_nu,g_aif,Npoints,N_ener,0);
		//cudaDeviceSynchronize();
		//Timer::getInstance().EndTimer("Kernal Execution");	
}

void GpuManager::ExecuteDopplerCrossSection(int N, int N_ener,int start_idx){
		cudaSetDevice(gpu_id);
		execute_two_step_kernal_doppler_block(g_freq, g_intens,g_energies,g_nu,g_gns,g_aif,N,N_ener,start_idx);
					//device_compute_cross_section_abscoef<<<gridSize,blockSize>>>(g_freq, g_cs,g_energies,g_gns,g_nu,g_aif,Npoints,N_ener,0);
		//cudaDeviceSynchronize();
		//Timer::getInstance().EndTimer("Kernal Execution");	
}

void GpuManager::ExecuteVoigtCrossSectionBlock(int N, int N_ener,int start_idx){
		cudaSetDevice(gpu_id);
		if(g_gamma==NULL)
			execute_two_step_kernal_voigt_block(g_freq, g_intens,g_energies,g_nu,g_gns,g_aif,N,N_ener,start_idx);
		else
			execute_two_step_kernal_voigt_block(g_freq, g_intens,g_energies,g_nu,g_gns,g_aif,g_gamma,g_n,N,N_ener,start_idx);
					//device_compute_cross_section_abscoef<<<gridSize,blockSize>>>(g_freq, g_cs,g_energies,g_gns,g_nu,g_aif,Npoints,N_ener,0);
		//cudaDeviceSynchronize();
		//Timer::getInstance().EndTimer("Kernal Execution");	
}


void GpuManager::ExecuteVoigtCrossSection(int N, int N_ener,int start_idx){
		cudaSetDevice(gpu_id);		
		execute_two_step_kernal_voigt(g_freq, g_intens,g_energies,g_nu,g_gns,g_aif,N,N_ener,start_idx);
					//device_compute_cross_section_abscoef<<<gridSize,blockSize>>>(g_freq, g_cs,g_energies,g_gns,g_nu,g_aif,Npoints,N_ener,0);
		//cudaDeviceSynchronize();
		//Timer::getInstance().EndTimer("Kernal Execution");	
}


void GpuManager::TransferResults(double* h_freq,double* h_intens,int N){
	cudaSetDevice(gpu_id);
	cudaDeviceSynchronize();
	cudaMemcpy(h_freq,g_freq,sizeof(double)*size_t(N),cudaMemcpyDeviceToHost);
	cudaMemcpy(h_intens,g_intens,sizeof(double)*size_t(N),cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	CheckCudaError("Copy Results");
}


void GpuManager::TrackMemory(size_t bytes){
	
	available_memory -= bytes;
}
void GpuManager::FreeMemory(size_t bytes){
	available_memory += bytes;
}
size_t GpuManager::GetAvailableMemory(){
	return available_memory;
}

void GpuManager::Cleanup(){
	cudaSetDevice(gpu_id);
	cudaFree(g_freq);
	cudaFree(g_intens);
	cudaFree(g_energies);
	cudaFree(g_gns);
	cudaFree(g_aif);
	cudaFree(g_gamma);
	cudaFree(g_n);


}	


