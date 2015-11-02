//#include "cuda_utils.cuh"
#include "Input.h"
#include "VoigtProfile.h"
#include "DopplerProfile.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>




int main(){
	printf("Hi!!\n");
	Timer::getInstance().StartTimer("Total Time");
	Input main_input;
	main_input.ReadInputII();
	BaseProfile* profile;	
	if(main_input.profile==VOIGT){
		profile = new VoigtProfile(&main_input);
	}else if (main_input.profile==DOPPLER){
		profile = new DopplerProfile(&main_input);
		
	}else{
		printf("Selected unimplemented profile!!\n");
		exit(0);
	}

	profile->Initialize();
	

	profile->ExecuteCrossSection();
	
	profile->OutputProfile();
	
	Timer::getInstance().EndTimer("Total Time");
	Timer::getInstance().PrintTimerInfo();

	//exit(0);
}
