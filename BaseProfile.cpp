#include "BaseProfile.h"

BaseProfile::BaseProfile(Input* pInput): input(pInput), Ntrans(0), h_energies(NULL), h_nu(NULL),h_gns(NULL), h_gammaL(NULL),h_n(NULL),h_aif(NULL) {



	start_nu = input->GetNuStart();
	end_nu = input->GetNuEnd();
	Npoints = input->GetNpoints();
	dfreq = (end_nu-start_nu)/double(Npoints);
	printf("Start_nu = %12.6f End_nu = %12.6f Dfreq = %12.6f Npoints =%d\n",start_nu,end_nu,dfreq,Npoints);

	freq = new double[Npoints];
	intens= new double[Npoints];


	for(int i = 0; i < Npoints; i++)
	{
		freq[i]=start_nu+double(i)*dfreq;
		intens[i]=0.0;
	}


}

BaseProfile::~BaseProfile(){

	delete[] freq;
	delete[] intens;

	


}

void BaseProfile::Initialize(){



	InitializeProfile();	
	//Get our transition files
	transition_files = input->GetTransFiles();
	
	

}


void BaseProfile::PerformCrossSection(double HW){
		printf("Waiting for previous calculation to finish......\n");
		fflush(0);
		Timer::getInstance().StartTimer("Execute GPU");	
		cudaDeviceSynchronize(); // Synchronize at the beginning
		gpu_manager->TransferVectors(Ntrans,h_energies, h_nu, h_aif,h_gns,h_gammaL,h_n);
		ib = std::max(round( ( min_nu-HW-start_nu)/dfreq ),0.0);
		ie =  std::min(round( ( max_nu+HW-start_nu)/dfreq ),double(Npoints));
		points = ie - ib;
		printf("Min nu = %12.6f Max_nu = %12.6f ib = %d ie = %d points = %d\n",min_nu,max_nu,ib,ie,points);
		fflush(0);
		gpu_manager->ExecuteCrossSection(points, Ntrans,ib);
		Timer::getInstance().EndTimer("Execute GPU");
}


void BaseProfile::OutputProfile(){
	//Output Data
	gpu_manager->TransferResults(freq,intens,Npoints);
	
	for(int i = 0; i < Npoints; i++){
		printf("%12.6f %13.8E\n",freq[i],intens[i]);
	}


}
