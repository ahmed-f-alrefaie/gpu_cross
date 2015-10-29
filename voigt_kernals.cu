//////////////////////////////////////////////////////////
/*			VOIGT PROFILE			*/
//////////////////////////////////////////////////////////
__global__ void device_compute_cross_section_voigt_stepone(double* g_energies,const int*  __restrict__ g_gns,const double*  __restrict__ g_nu,const double*  __restrict__ g_aif, const int N_ener){
	//The stored shared data
	
	//Get the global and local thread number
	int g_idx = blockIdx.x * blockDim.x + threadIdx.x;
	double ei,gns,nu_if,aif,abscoef;
	double temp_2 = 1.0;//cross_constants.ln2pi/cross_constants.halfwidth;
	
	//if(g_idx == 0) printf("partition = %12.6f\n",cross_constants.partition);
	if(g_idx < N_ener){
			//Store values in local memory
			ei = g_energies[g_idx];
			gns = g_gns[g_idx];
			nu_if = g_nu[g_idx];
			aif = g_aif[g_idx];
				
			abscoef= cross_constants.cmcoef*temp_2*aif*gns
				*exp(-cross_constants.beta*ei)*(1.0-exp(-cross_constants.beta*nu_if))/
				(nu_if*nu_if*cross_constants.partition);
			if(nu_if==0)abscoef=0.0;
			g_energies[g_idx] = abscoef;
			
	}


}

__global__ void device_compute_cross_section_voigt_stepone(double* g_energies,const int*  g_gns,const double*  g_nu,const double*  g_aif,double*  g_gamma,double*  g_n, const int N_ener){
	//The stored shared data
	
	//Get the global and local thread number
	int g_idx = blockIdx.x * blockDim.x + threadIdx.x;
	double ei,gns,nu_if,aif,abscoef;
	double gammaL;
	
	//cross_constants.ln2pi/cross_constants.halfwidth;
	//if(g_idx == 0) printf("partition = %12.6f\n",cross_constants.partition);
	if(g_idx < N_ener){
			//Store values in local memory
			ei = g_energies[g_idx];
			gns = g_gns[g_idx];
			nu_if = g_nu[g_idx];
			aif = g_aif[g_idx];

			if(nu_if==0) nu_if = 1e-6;
			abscoef= cross_constants.cmcoef*aif*gns
				*exp(-cross_constants.beta*ei)*(1.0-exp(-cross_constants.beta*nu_if))/
				(nu_if*nu_if*cross_constants.partition);
			if(gns==-1) abscoef = aif;
			g_energies[g_idx] = abscoef;

			gammaL = g_gamma[g_idx]*pow(296.0/cross_constants.temperature,g_n[g_idx])*cross_constants.pressure; 
			g_gamma[g_idx] = gammaL;

			//if(threadIdx.x == 0) printf("%14.2E   %14.2E\n",abscoef,gammaL) ;
			
	}


}


__global__ void device_compute_cross_section_voigt_steptwo(const double* g_freq, double* g_cs,const double*  g_nu,const double*  g_abscoef,const int N,const int N_ener,const int start_idx){
	//The stored shared data
	__shared__ double l_nu[VOIGT_SHARED_SIZE];
	__shared__ double l_abscoef[VOIGT_SHARED_SIZE];
	//Get the global and local thread number
	int g_idx = blockIdx.x * blockDim.x + threadIdx.x;
	int l_idx = threadIdx.x;
	int block_dim = VOIGT_SHARED_SIZE;
	double cs_val = 0.0;
	double dfreq_=0.0;
	double freq = 0.0;
	double gammaG=0.05,gammaL=0.05,x,y;
	double dpwcoeff = sqrt(2.0*LN2*BOLTZ*cross_constants.temperature/(cross_constants.mean_mass))/VELLGT;
	//double nu_if;
	//if(g_idx == 0) printf("BLOCK_SIZE = %d\n",blockDim.x);
	if(g_idx < N){
		freq = g_freq[start_idx+g_idx];
		//cs_val = g_cs[start_idx+g_idx];
	}

	//if(g_idx==9999)  printf("%12.6f\n",freq);	

	for(int i = 0; i < N_ener; i+=VOIGT_SHARED_SIZE){
		l_nu[l_idx] = 1.0;
		l_abscoef[l_idx] = 0.0;

		if(i + l_idx < N_ener)
		{	
			l_nu[l_idx] = g_nu[i + l_idx];
			l_abscoef[l_idx] = g_abscoef[i + l_idx];
		}
		
		__syncthreads();
		for(int j = 0; j < VOIGT_SHARED_SIZE; j++){
			dfreq_=l_nu[j]-freq;
			gammaG = l_nu[j]*dpwcoeff;
			x =SQRTLN2*abs(dfreq_)/gammaG;
			y =SQRTLN2*gammaL/gammaG;
			double xxyy = x * x + y * y;

			//Algorithm 916
			if(xxyy < 100.0){
				cs_val+=l_abscoef[j]*SQRTLN2PI/(gammaG)*y*voigt_916(x,y,1.0);					
			}
else{
				//3-point gauss hermite
			cs_val+=l_abscoef[j]*(SQRTLN2PI/gammaG)*voigt_threegausshermite(x,y,xxyy);
			}
			//*__expf(temp_3*dfreq_*dfreq_);

			
		}
		__syncthreads();
		


	}
	

	if(g_idx < N) g_cs[start_idx+g_idx]+=cs_val;


}




__global__ void device_compute_cross_section_voigt_steptwo_block(const double*  g_freq, double* g_cs,const double*   g_nu,const double*  g_abscoef,const int N,const int N_ener,const int start_idx){
	//The stored shared data
	//__shared__ double l_nu[BLOCK_SIZE];
	//__shared__ double l_abscoef[BLOCK_SIZE];
	__shared__ double l_cs_result[VOIGT_BLOCK];
	//Get the global and local thread number
	int b_idx = blockIdx.x;
	int l_idx = threadIdx.x;
	double cs_val = 0.0;
	double dfreq_=0.0;
	double freq = 0.0;
	double nu = 0.0;
	double gammaG=0.05,gammaL=0.05;
	double x,y;

	double dpwcoeff = sqrt(2.0*LN2*BOLTZ*cross_constants.temperature/(cross_constants.mean_mass))/VELLGT;

	//double temp_2=cross_constants.ln2pi/cross_constants.halfwidth;
	//double temp_3 = -cross_constants.ln2*(1.0/(cross_constants.halfwidth*cross_constants.halfwidth));

	freq = g_freq[start_idx + b_idx];
	//cs_val = g_cs[start_idx+g_idx];

	//if(g_idx==9999)  printf("%12.6f\n",freq);	
	l_cs_result[l_idx] = cs_val;
	for(int i = l_idx; i < N_ener; i+=VOIGT_BLOCK){
		nu = 0.0;
		//Read value of nu
		nu = g_nu[i];
		dfreq_ = nu-freq;
		if(dfreq_ < -500.0*gammaL)
			continue;
		if(dfreq_ > 500.0*gammaL)
			break;
		gammaG = nu*dpwcoeff;
		x =SQRTLN2*dfreq_/gammaG;
		y =SQRTLN2*gammaL/gammaG;
		double xxyy = x * x + y * y;

		
		////Algorithm 916
		if(xxyy < 100.0){
			cs_val+=g_abscoef[i]*SQRTLN2PI/(gammaG)*y*voigt_916(x,y,1.0);					
		}else{
			//3-point gauss hermite
			cs_val+=g_abscoef[i]*(SQRTLN2PI/gammaG)*voigt_threegausshermite(x,y,xxyy);
		}
			

	}
	//Store results into shared memory
	l_cs_result[l_idx] = cs_val;
	cs_val = 0;
	//Wait for everone to finish nicely
	__syncthreads();
	if(l_idx == 0){
		for(int i = 0; i < VOIGT_BLOCK; i++)
			cs_val+=l_cs_result[i];
		
		g_cs[start_idx+b_idx]+=cs_val;		
	}

}

__global__ void device_compute_cross_section_voigt_steptwo_block(const double*  g_freq, double* g_cs,const double*   g_nu,const double*  g_abscoef,const double*  g_gamma,const int N,const int N_ener,const int start_idx){
	//The stored shared data
	//__shared__ double l_nu[BLOCK_SIZE];
	//__shared__ double l_abscoef[BLOCK_SIZE];
	__shared__ double l_cs_result[VOIGT_BLOCK];
	//Get the global and local thread number
	int b_idx = blockIdx.x;
	int l_idx = threadIdx.x;
	double cs_val = 0.0;
	double dfreq_=0.0;
	double freq = 0.0;
	double nu = 0.0;
	double gammaG=0.05,gammaL=0.05;
	double x,y;

	double dpwcoeff = sqrt(2.0*BOLTZ*cross_constants.temperature*NA/((cross_constants.mean_mass)))/VELLGT;

	//double temp_2=cross_constants.ln2pi/cross_constants.halfwidth;
	//double temp_3 = -cross_constants.ln2*(1.0/(cross_constants.halfwidth*cross_constants.halfwidth));

	freq = g_freq[start_idx + b_idx];
	//cs_val = g_cs[start_idx+g_idx];

	//if(g_idx==9999)  printf("%12.6f\n",freq);	
	l_cs_result[l_idx] = cs_val;
	for(int i = l_idx; i < N_ener; i+=VOIGT_BLOCK){
		nu = 0.0;
		//Read value of nu
		nu = g_nu[i];
		dfreq_ = nu-freq;
		gammaL = g_gamma[i];
		if(dfreq_ < -500.0*gammaL)
			continue;
		if(dfreq_ > 500.0*gammaL)
			break;
		gammaG = 1.0/(nu*dpwcoeff);
		x =abs(dfreq_)*gammaG;
		y =gammaL*gammaG;
		double xxyy = x * x + y * y;
		double voigt_check;// = voigt_916(x,y,0.9);


		
		////Algorithm 916
		if(xxyy < 100.0){
			voigt_check = voigt_916(x,y,0.9);
			//cs_val+=g_abscoef[i]*voigt_check*gammaG*ISQRTPI;					
		}else if(xxyy < 1.0e6){
			//3-point gauss hermite
			voigt_check = voigt_threegausshermite(x,y,xxyy);
			//cs_val+=g_abscoef[i]*ISQRTPI*gammaG;
		}else{
			voigt_check = y/(PI*xxyy);
			//cs_val+= g_abscoef[i]*ISQRTPI*gammaG;
		}
		cs_val+=g_abscoef[i]*voigt_check*gammaG*ISQRTPI;
		//if((blockIdx.x * blockDim.x + threadIdx.x)==0)  printf("dfreq = %14.4E x=%14.4E y=%14.4E gammaL = %14.4E gammaG = %14.4E abscoef=%14.4E voigt=%14.4E cs_val=%14.4E\n",dfreq_,x,y,gammaL,gammaG,g_abscoef[i],voigt_check,cs_val);			

	}
	//Store results into shared memory
	l_cs_result[l_idx] = cs_val;
	cs_val = 0;
	//Wait for everone to finish nicely
	__syncthreads();
	if(l_idx == 0){
		for(int i = 0; i < VOIGT_BLOCK; i++)
			cs_val+=l_cs_result[i];
		
		g_cs[start_idx+b_idx]+=cs_val;		
	}

}
