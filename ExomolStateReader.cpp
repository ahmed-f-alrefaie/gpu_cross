#include "ExomolStateReader.h"
#include <fstream>
#include <sstream>
#include <map>
using namespace std;

bool ExomolStateReader::GetGammaN(state* state_i,state* state_f, double & gamma, double & n){
	if(broad_info.broad_map.size()==0)
		return false;
	int Ji = state_i->J;
	int Jf = state_f->J;
	std::stringstream key;
	key<<"a1"<<Ji<<Jf;
	if(broad_info.broad_map.count(key.str()) > 0){
		gamma = broad_info.broad_map[key.str()].gamma;
		n = broad_info.broad_map[key.str()].n;
		gamma*= m_pressure;
		return true;
		
	}
	key.clear();
	key<<"a0"<<Ji<<-1;
	if(broad_info.broad_map.count(key.str())>0){
		gamma = broad_info.broad_map[key.str()].gamma;
		n = broad_info.broad_map[key.str()].n;
		gamma*= m_pressure;
		return true;
		
	}

	return false;
		
	



}

void ExomolStateReader::InitializeBroadener(std::string broadener){
	std::ifstream i_broad(broadener.c_str());
	if(!i_broad){
		std::cout<<"Broadener file with name "<<broadener<<" not found!"<<std::endl;
		use_broadeners = false;
		return;
	}
	
	double gam;
	double n;
	std::string line;
	int Ji,Jf;
	string code;
	while(getline(i_broad,line)){
		//
		code = "a0";
		vector<string> split_line = split(line);
		Ji=atoi(trim(split_line[3]).c_str());
		Jf=-1;
		gam=atof(trim(split_line[1]).c_str());
		n=atof(trim(split_line[2]).c_str());
		if(trim(split_line[0])=="a1"){
			code="a1";
			Jf=atoi(trim(split_line[4]).c_str());	
		}
		//Construct the code
		std::stringstream key;
		key<<code<<Ji<<Jf;
		std::cout<<"Key: "<<key.str()<<std::endl;
		//broadData broad;
		broad_info.broad_map[key.str()].gamma = gam;
		broad_info.broad_map[key.str()].n = n;
				
		//broad_info.broad_map[key.str()]= broad;
		
		
		
		


	}	
	use_broadeners = true;
	//Begin reading
	




}
ExomolStateReader::ExomolStateReader(std::string pFilename,double partition,double pressure,std::string pbroadener) : StateReader(pFilename,partition) ,m_pressure(pressure){

	printf("Using Exomol File format\n");
	initialize_states(pFilename.c_str(), &exomol);
	
	
	//Initializwe the broadeners
	InitializeBroadener(pbroadener);
	


}

ExomolStateReader::~ExomolStateReader(){
	CloseFile();
}

bool ExomolStateReader::OpenFile(std::string pFilename){
	stream.open(pFilename.c_str());
}
bool ExomolStateReader::CloseFile(){
	stream.close();
}
bool ExomolStateReader::ReadNextState(double & nu,double & gns,double & e_i, double & aif, double & gam,double & n){
	int id_f,id_i;	
	if ( (stream>>id_f>>id_i>>aif)==false)
		return false;
	
	id_f-=1;
	id_i-=1;

	nu = exomol.states[id_f].energy-exomol.states[id_i].energy;
	e_i = exomol.states[id_i].energy;
	gns = exomol.states[id_f].gns;
						//if(aif < 1e-30) continue;
	if(use_broadeners=true)
		GetGammaN(&exomol.states[id_i],&exomol.states[id_f],gam,n);
	
	//printf("%12.6f %12.6f %d %14.3E\n",nu,e_i,gns,gam,n);
	return true;	

}
bool ExomolStateReader::ReadNextState(double & nu,double & gns,double & e_i, double & aif){
	double dum_g=0,dum_n=0;
	ReadNextState(nu,gns,e_i, aif, dum_g,dum_n);
}
double ExomolStateReader::ComputePartition(double temperature){
	if(partition<0.0){
		this->partition = compute_partition(&exomol,temperature);
		std::cout<<"Partition is: "<<partition<<std::endl;
	}
	return partition;
}

