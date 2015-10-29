#include "MultiGPUManager.h"


MultiGpuManager::MultiGpuManager(ProfileType pProfile) : BaseManager(pProfile), selected_gpu(0){
	int devCount;
	cudaGetDeviceCount(&devCount);
	printf("CUDA Device Query...\n");
	printf("There are %d CUDA devices.\n", devCount);

	t_num_gpus=devCount;

    	//Alocated all the managers needed to do gpu shindig
	gpu_managers = new GpuManager*[t_num_gpus];
	for(int i = 0; i < t_num_gpus; i++){
		gpu_managers[i] = new GpuManager(pProfile,i);
	}
	
	selected_manager = gpu_managers[selected_gpu];
}




MultiGpuManager::~MultiGpuManager(){
	Cleanup();
}
void MultiGpuManager::InitializeVectors(int Npoints){
	N_trans = 1000000000;
	for(int i = 0; i < t_num_gpus; i++){
		gpu_managers[i]->InitializeVectors(Npoints);
		// get the largest NTrans
	}
		
	N_trans = gpu_managers[0]->GetNtrans(); 
	printf("MultiGPUManger: N_trans = %d\n",N_trans);

}
void MultiGpuManager::InitializeConstants(double half_width,double temperature, double partition,double dfreq,double meanmass,double pressure,double ref_temp){

	for(int i = 0; i < t_num_gpus; i++){
		gpu_managers[i]->InitializeConstants(half_width,temperature, partition,dfreq,meanmass,pressure,ref_temp);
	}


}
void MultiGpuManager::TransferVectors(size_t Nener,double* h_energies, double* h_nu, double* h_aif,int* h_gns,double* h_gamma,double * h_n){

	selected_manager->TransferVectors(Nener,h_energies, h_nu, h_aif,h_gns,h_gamma,h_n);

}
void MultiGpuManager::TransferFreq(double* h_freq,double* h_intens,int N){
	for(int i = 0; i < t_num_gpus; i++){
		gpu_managers[i]->TransferFreq( h_freq,h_intens,N);
	}	

}
void MultiGpuManager::ExecuteCrossSection(int N, int N_ener,int start_idx){
	printf("Executing cross-section on GPU %d\n",selected_gpu);	
	selected_manager->ExecuteCrossSection(N,N_ener,start_idx);
	SwitchGPU(); //Switch to the next gpu
}
void MultiGpuManager::TransferResults(double* h_freq,double* h_intens,int N){
	double* temp_freq = new double[N];
	double* temp_intens = new double[N];
	for(int i = 0; i < N; i++){
		temp_intens[i] = 0.0;
	}	
	printf("Transferring results from all GPUs into host!\n");
	
	for(int i = 0; i < t_num_gpus; i++){
		gpu_managers[i]->TransferResults(temp_freq,temp_intens,N);
		for(int i = 0; i < N; i++){
			h_intens[i] += temp_intens[i];
		}			


	}

	delete[] temp_freq;
	delete[] temp_intens;
	


}
void MultiGpuManager::Cleanup(){
	for(int i = 0; i < t_num_gpus; i++){
		gpu_managers[i]->Cleanup();
	}	
	delete[] gpu_managers;

}
void MultiGpuManager::SwitchGPU(){
	selected_gpu++;
	if(selected_gpu>=t_num_gpus)
		selected_gpu = 0;
	
	selected_manager = gpu_managers[selected_gpu];
}


