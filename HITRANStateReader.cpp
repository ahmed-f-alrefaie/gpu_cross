#include "HITRANStateReader.h"

int GetFieldStart(int field){
	const int field_lengths[] =  {2,1,12,10,10,5,5,10,4,8,15,15,15,15,6,12,1,7,7};
	int start = 0;	
	for(int i = 0; i < field; i++){
		start=start+field_lengths[i];
	} 
		
	return start;
	
}

int GetFieldLength(int field){
	const int field_lengths[] =  {2,1,12,10,10,5,5,10,4,8,15,15,15,15,6,12,1,7,7};
	return field_lengths[field];
	
}


HITRANStateReader::HITRANStateReader(std::string pFilename,double pPartition,double pressure,double temperature,double mixture):StateReader(pFilename,pPartition) , m_pressure(pressure),m_mixture(mixture), m_temperature(temperature){
	 open_file = OpenFile("");

}

bool HITRANStateReader::OpenFile(std::string pFilename){
	if(pFilename == "")
		return false;
	if(open_file==true){
		CloseFile();
	}
	
	hitran_file = fopen(m_filename.c_str(),"r");
	if(hitran_file==NULL){
		printf("Could not open file %s ",m_filename.c_str());
		return false;
	}
	open_file = true;
	return true;
	
}

bool HITRANStateReader::CloseFile(){
	if(open_file==true && hitran_file != NULL){
		fclose(hitran_file);
	}


}

bool HITRANStateReader::ReadNextState(double & nu,double & gns,double & e_i, double & aif){

	double dummy_gam,dummy_n;
	return ReadNextState(nu,gns,e_i,aif,dummy_gam,dummy_n);

	
	





}
bool HITRANStateReader::ReadNextState(double & nu,double & gns,double & e_i, double & aif, double & gam,double & n){

	if(!open_file)
		return false;

	char c_nu[13];
	char c_eif[13];
	char c_gns[13];
	char c_aif[13];
	char c_gam[13];
	char c_n[13];
	char c_self[13];
	
	//printf("%s\n",buffer);
	if(NULL==fgets(buffer,1024,hitran_file))
		return false;
	
	//for(int i = 0; i < 160; i++)
	//	if(buffer[i]==' ') buffer[i]='0';

//	scanf(buffer,"%2u%1u%12f%10f%10f%5f%5f%10f%4f%8f%15c%15c%15c%15c%6c%12c%1c%7f%7f",
//&d1,&d1,&nu,&dd1,&aif,&dd1,&dd1,&e_i,&dd1,&dd1,small_buff,small_buff,small_buff,
//	small_buff,small_buff,&d1,&d1,small_buff,&gns,&dd1);
	int start = GetFieldStart(2);
	int length = GetFieldLength(2);


	memcpy(c_nu,buffer+start,length);
	start = GetFieldStart(4);
	length = GetFieldLength(4);
	memcpy(c_aif,buffer+start,length);	
		
	start = GetFieldStart(7);
	length = GetFieldLength(7);
	memcpy(c_eif,buffer+start,length);	

	start = GetFieldStart(5);
	length = GetFieldLength(5);
	memcpy(c_gam,buffer+start,length);	

	start = GetFieldStart(6);
	length = GetFieldLength(6);
	memcpy(c_self,buffer+start,length);	

	start = GetFieldStart(8);
	length = GetFieldLength(8);
	memcpy(c_n,buffer+start,length);	

	start = GetFieldStart(17);
	length = GetFieldLength(17);
	memcpy(c_gns,buffer+start,length);

	

	c_nu[12]='\0';
	c_eif[12]='\0';
	c_aif[12]='\0';
	c_gns[12]='\0';
	c_gam[12]='\0';
	c_n[12] = '\0';
	c_self[12] = '\0';
	//printf("Here: %s %14.3E\n",c_nu,atof(c_nu));
	nu = atof(c_nu);
	e_i = atof(c_eif);
	aif=atof(c_aif);
	gns = atof(c_gns);
	gam = (atof(c_gam)*m_mixture + atof(c_self)*(1.0-m_mixture))*m_pressure;
	n=atof(c_n);

	if(aif==0.0) {
		//Compute the 'GNSAIF'
		
		start = GetFieldStart(3);
		length = GetFieldLength(3);
		memcpy(c_aif,buffer+start,length);	
		c_aif[12]='\0';
		//printf("%s\n",buffer);
		//printf("%12.6f %12.6f %14.4E %d %12.6f %12.6f\n",nu, e_i, aif, int(gns),gam,n);
		aif = atof(c_aif);
		gns=-1.0;
		
	}
	return true;

}

HITRANStateReader::~HITRANStateReader(){
	CloseFile();
}
