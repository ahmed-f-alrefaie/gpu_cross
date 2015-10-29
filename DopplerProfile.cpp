#include "DopplerProfile.h"
#include "Timer.h"

const double dopplerHW = 100.0;

void DopplerProfile::InitializeProfile(){
	
	printf("Profile Selected is Doppler\n");
	//Initialize the stateReader TODO: Use Factory Design Pattern for this
	if(input->GetFileType()==HITRAN_TYPE){
		state_reader = (StateReader*)(new HITRANStateReader( input->GetTransFiles()[0] , input->GetPartition(),input->GetPressure(),input->GetTemperature(),input->GetHitranMixture()   ));
	}else if(input->GetFileType()==EXOMOL_TYPE){
		state_reader = (StateReader*)(new ExomolStateReader( input->GetStateFile() , input->GetPartition(),input->GetPressure()));
	//	exit(0);
	}else{
		exit(0);
	}

	gpu_manager=new MultiGpuManager(profile);
	gpu_manager->InitializeVectors(Npoints);
	gpu_manager->TransferFreq(freq,intens,Npoints);
	num_trans_fit=gpu_manager->GetNtrans();


	gpu_manager->InitializeConstants(input->GetHalfWidth(),input->GetTemperature(), state_reader->ComputePartition(input->GetTemperature()),dfreq,input->GetMeanMass(),input->GetPressure(),296.0);


	h_energies = new double[num_trans_fit];
	h_nu = new double[num_trans_fit];
	h_aif = new double[num_trans_fit];
	h_gns = new int[num_trans_fit];
	
}

void DopplerProfile::ExecuteCrossSection(){
	//We have no transitions right now	
	Ntrans = 0;
	double t_nu,t_ei,t_aif,t_gns;

	int num_trans_files = transition_files.size();
	max_nu = 0.0;
	min_nu = 99999.0;

	printf("Range = %12.6f - %12.6f\n",start_nu-dopplerHW,end_nu+dopplerHW);
	for(int  i = 0; i < num_trans_files; i++){
		state_reader->OpenFile(transition_files[i]);
		printf("%s\n",transition_files[i].c_str());
		fflush(0);
		Timer::getInstance().StartTimer("FileIO");
		while(state_reader->ReadNextState(t_nu,t_gns,t_ei,t_aif)){


			
			if(t_nu< start_nu-dopplerHW)
				continue;
			if(t_nu> end_nu+dopplerHW)
				break;

			h_energies[Ntrans] = t_ei;
			h_nu[Ntrans] = t_nu;
			h_gns[Ntrans] = int(t_gns);
			h_aif[Ntrans] = t_aif;				
			
			
			//Find the maximum and minimum nu 
			min_nu=std::min(min_nu,h_nu[Ntrans]);
			max_nu=std::max(max_nu,h_nu[Ntrans]);		
			//max_gammaL = std::max(max_gammaL,h_gammaL[Ntrans]);
			//min_n = std::min(min_n,h_n[Ntrans]); 
			//Restrict to a set number of points
			
			//if(
			Ntrans++;
			if(Ntrans>=num_trans_fit){
				//printf("%d\n",temp_points);
				//gammaHW = 500.0*max_gammaL*pow((300.0/input->GetTemperature()),min_n)*input->GetPressure(); 
				PerformCrossSection(dopplerHW);


				//Reset counters		
				max_nu = 0.0;
				min_nu = 99999.0;
				Ntrans = 0;
			}
		}
		Timer::getInstance().EndTimer("FileIO");
		state_reader->CloseFile();

	}
	printf("Left over transitions: %d\n",Ntrans);
	if(Ntrans > 0 ){

		//gammaHW = 500.0*max_gammaL*pow((300.0/input->GetTemperature()),min_n)*input->GetPressure(); 
		PerformCrossSection(dopplerHW);

				//Reset counters		
	}



}

