#include "BaseProfile.h"
#include "HITRANStateReader.h"
#include "ExomolStateReader.h"
#include <cmath>
#pragma once

class DopplerProfile : public BaseProfile{	
	private:
		//The double and	
		double DopplerCoeff;

	public:
		DopplerProfile(Input* pInput) : BaseProfile(pInput){	profile = DOPPLER;};
		void InitializeProfile();
		void ExecuteCrossSection();


};


