#include "OpenMPManager.h"
#include "defines.h"
#include <cstring>
#include "profiles.h"
		
OpenMPManager::OpenMPManager(ProfileType pprofile, int num_threads,size_t total_memory): BaseManager(pprofile), total_threads(num_threads){
	printf("Utilizing Threading\n");
	//Initialize memory
	InitializeMemory(total_memory/2l); //Half the memory we specify
	//Initialize our threads
	//workers = new std::thread[total_threads];
	t_intens = new double*[total_threads];
	threads_done = new bool[total_threads];
	for(int i = 0; i < total_threads; i++)
		threads_done[i] = true;
	

}
void OpenMPManager::InitializeVectors(int Npoints){
	g_freq = new double[Npoints];
	TrackMemory(size_t(Npoints)*sizeof(double));

	g_intens = new double[Npoints];
	TrackMemory(size_t(Npoints)*sizeof(double));
	
	//We create seperate intensities for each thread
	for(int i = 0; i < total_threads; i++){
		t_intens[i] = new double[Npoints];
		TrackMemory(size_t(Npoints)*sizeof(double));
	}

	if(profile==VOIGT)
		N_trans = GetAvailableMemory()/(sizeof(double)+sizeof(double)+sizeof(double)+sizeof(int)+sizeof(double) + sizeof(double));
	else
		N_trans = GetAvailableMemory()/(sizeof(double)+sizeof(double)+sizeof(double)+sizeof(int));

	printf("Number of transitions in memory = %d\n",N_trans);

	g_energies = new double[N_trans];
	TrackMemory(sizeof(double)*size_t(N_trans));
	g_nu = new double[N_trans];
	TrackMemory(sizeof(double)*size_t(N_trans));
	g_aif = new double[N_trans];
	TrackMemory(sizeof(double)*size_t(N_trans));
	g_gns = new int[N_trans];	
	TrackMemory(sizeof(int)*size_t(N_trans));

	if(profile==VOIGT){
		g_gamma = new double[N_trans];
		TrackMemory(sizeof(double)*size_t(N_trans));
		g_n = new double[N_trans];
		TrackMemory(sizeof(double)*size_t(N_trans));

	}

}
void OpenMPManager::InitializeConstants(double half_width,double temperature, double partition,double dfreq,double meanmass,double pressure,double ref_temp){
		m_half_width=half_width;
		m_temperature=temperature;
		m_partition=partition;
		m_dfreq=dfreq;
		m_mean_mass=meanmass;
		m_pressure=pressure;
		m_ref_temp=ref_temp;
		m_cmcoef = 1.0/(8.0*PI*VELLGT); 
		m_beta = PLANCK*VELLGT/(BOLTZ*temperature);
}
void OpenMPManager::TransferVectors(size_t Nener,double* h_energies, double* h_nu, double* h_aif,int* h_gns,double* h_gamma,double * h_n){

	JoinAllThreads();
	memcpy(g_energies,h_energies,size_t(Nener)*sizeof(double));
	memcpy(g_nu,h_nu,size_t(Nener)*sizeof(double));
	memcpy(g_aif,h_aif,size_t(Nener)*sizeof(double));
	memcpy(g_gns,h_gns,size_t(Nener)*sizeof(int));
	if(h_gamma!= NULL && h_n != NULL){
		memcpy(g_gamma,h_gamma,size_t(Nener)*sizeof(double));
		memcpy(g_n,h_n,size_t(Nener)*sizeof(double));
	} 
	
}

void OpenMPManager::TransferFreq(double* h_freq,double* h_intens,int N){
	memcpy(g_freq,h_freq,size_t(N)*sizeof(double));
	for(int i = 0 ; i < total_threads; i++){
		memcpy(t_intens[i],h_intens,sizeof(double)*size_t(N));
	}
	start_nu = g_freq[0];
}
void OpenMPManager::ExecuteCrossSection(int N, int N_ener,int start_idx){
	printf("Executing cross-section with %d threads\n",total_threads);
	//printf("Here we go!\n");
	///fflush(0);
	if(profile == VOIGT){
		//if(N > 100000)
		//	ExecuteVoigtCrossSection(N, N_ener,start_idx);
		//else
		for(int i = 0; i < total_threads; i++){
			// printf("PUSH");
			//fflush(0);
			 threads_done[i]=false;
			// printf("HARD");	
			//fflush(0);		
			 workers.push_back(std::thread(&OpenMPManager::ComputeVoigt,this,i,N, N_ener,start_idx));
		}
	}else if(profile == DOPPLER){
		for(int i = 0; i < total_threads; i++){
			 threads_done[i]=false;		
			 workers.push_back(std::thread(&OpenMPManager::ComputeDoppler,this,i,N, N_ener,start_idx));
		}
	}else{
		printf("Profile not implemented!!\n");
		exit(0);
	}
	//printf("Started_threads\n");
	fflush(0);
	
}
void OpenMPManager::TransferResults(double* h_freq,double* h_intens,int N){

	printf("Wait for the jobs to finish\n");
	JoinAllThreads(); // If we aren't done wait
	printf("Transferring from all threads into host!\n");
	for(int i = 0; i < total_threads; i++){
		for(int j = 0; j < N; j++){
			g_intens[j]+=t_intens[i][j];
		}
	}

	memcpy(h_intens,g_intens,sizeof(double)*size_t(N));


}
void OpenMPManager::Cleanup()
{

}

void OpenMPManager::ComputeDoppler(int id,int Npoints,int Nener,int start_idx){	

	int thread_id = id;
	int for_block = std::ceil(float(Nener)/float(total_threads));
	int start = id*for_block;
	int end = std::min(id*for_block + for_block,Nener);
	//printf("Thread id = %d start = %d, end = %d block_size=%d\n",thread_id,start,end,for_block);
	double* intens_ptr = t_intens[thread_id];
	double dpwcoef = sqrt(2.0*LN2*BOLTZ*m_temperature*NA/(m_mean_mass))/VELLGT;
	double lorentz_cutoff = min(25.0*m_pressure,100.0);
	double x0 = m_dfreq*0.5;
	//fflush(0);
	for(int i = start; i < end; i++){
		double nu_if = g_nu[i];
		
		if(nu_if==0) nu_if = 1e-6;

		double abscoef= m_cmcoef*g_aif[i]*g_gns[i]
				*exp(-m_beta*g_energies[i])*(1.0-exp(-m_beta*nu_if))/
				(nu_if*nu_if*m_partition);	
		
		//double gammaL = g_gamma[i]*pow(m_ref_temp/m_temperature,g_n[i]);
		int ib = std::max(round( ( nu_if-DOPPLER_CUTOFF-start_nu)/m_dfreq ),double(start_idx));
		int ie =  std::min(round( ( nu_if+DOPPLER_CUTOFF-start_nu)/m_dfreq ),double(Npoints + start_idx));
		//printf("[%d] nu = %12.6f abscoef = %14.3E ib = %d, ie = %d\n",thread_id,nu_if,abscoef,ib,ie);
		if(abscoef==0.0) continue;
		for(int j = ib; j < ie; j++){
			double gammaG = 1.0/(nu_if*dpwcoef);
			double dfreq_ = abs(nu_if - g_freq[j]);
			double xp,xm,de;
			xp = gammaG*((dfreq_)+x0);
			xm = gammaG*((dfreq_)-x0);
			de = erf(xp)-erf(xm); 
			intens_ptr[j]+=abscoef*de;
			//printf("[%d] x=%14.3E y = %14.3E Humlik = %14.3E\n",thread_id,x,y,humlic(x,y));
		}

	
	}
	threads_done[thread_id]=true;	
	


}	

void OpenMPManager::ComputeVoigt(int id,int Npoints,int Nener,int start_idx){
	
	int thread_id = id;
	int for_block = std::ceil(float(Nener)/float(total_threads));
	int start = id*for_block;
	int end = std::min(id*for_block + for_block,Nener);
	//printf("Thread id = %d start = %d, end = %d block_size=%d\n",thread_id,start,end,for_block);
	double* intens_ptr = t_intens[thread_id];
	double dpwcoef = sqrt(2.0*LN2*BOLTZ*m_temperature*NA/(m_mean_mass))/VELLGT;
	double lorentz_cutoff = min(25.0*m_pressure,100.0);
	
	//fflush(0);
	for(int i = start; i < end; i++){
		double nu_if = g_nu[i];
		
		if(nu_if==0) nu_if = 1e-6;

		double abscoef= m_cmcoef*g_aif[i]*g_gns[i]*ISQRTPI
				*exp(-m_beta*g_energies[i])*(1.0-exp(-m_beta*nu_if))/
				(nu_if*nu_if*m_partition);	
		
		double gammaL = g_gamma[i]*pow(m_ref_temp/m_temperature,g_n[i]);
		int ib = std::max(round( ( nu_if-lorentz_cutoff-start_nu)/m_dfreq ),double(start_idx));
		int ie =  std::min(round( ( nu_if+lorentz_cutoff-start_nu)/m_dfreq ),double(Npoints + start_idx));
		//printf("[%d] nu = %12.6f abscoef = %14.3E ib = %d, ie = %d\n",thread_id,nu_if,abscoef,ib,ie);
		if(abscoef==0.0) continue;
		for(int j = ib; j < ie; j++){
			double gammaG = 1.0/(nu_if*dpwcoef);
			double dfreq_ = nu_if - g_freq[j];
			double x =abs(dfreq_)*gammaG;
			double y =gammaL*gammaG;
			intens_ptr[j]+=abscoef*humlic(x,y)*gammaG;
			//printf("[%d] x=%14.3E y = %14.3E Humlik = %14.3E\n",thread_id,x,y,humlic(x,y));
		}

	
	}
	threads_done[thread_id]=true;	
	


}

bool OpenMPManager::ReadyForWork(){

	bool status=true;
	for(int i = 0; i < total_threads; i++){
		status &= threads_done[i];
	}
	return status;
}

void OpenMPManager::JoinAllThreads(){
	printf("Joining threads...");
	for(int i = 0; i < workers.size();i++){
		//printf("..%d..",i);
		workers[i].join();
	}
	workers.clear(); // Destroy the threads
	printf("done!\n");
}
