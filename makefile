goal:   gpu_cross.x



PLAT = k20
FOR  = nvcc 

#FOR = /opt/intel/composer_xe_2011_sp1.11.339/bin/intel64/ifort
FFLAGS = --ptxas-options=-v -O3 -arch=sm_35 -Xptxas -v -lineinfo
#FFLAGS =  -O3 -xHost -align -ansi-alias -mcmodel=medium -g -traceback  -openmp
###FFLAGS = -O3 -openmp -g -traceback -xHost -align -ansi-alias

###PLASMA = -lplasma -lcoreblas -lquark -llapacke -lpthread -lhwloc -lmrrr -lhwloc  -lirc

####PLASMA =-lplasma -lcoreblasqw -lcoreblas -lquark      -llapacke -lpthread -lhwloc -lmrrr -lhwloc  -lirc

#PLIBS  =  ../../plasma/libplasma.a ../../plasma/libcoreblas.a ../../plasma/libquark.a ../../plasma/libmrrr.a 

#LIBS = -lplasma -lcoreblasqw -lcoreblas -lquark     -llapacke -mkl=parallel -lpthread -lhwloc -lmrrr -lhwloc  -lirc


FC = ifort
CC = g++
CCFLAGS = -O3
##FC=/sw/sdev/intel/composer_xe_2013.1.117/composer_xe_2013/bin/ifort
##CC=/sw/sdev/intel/composer_xe_2013.1.117/composer_xe_2013/bin/icc

#FC=/opt/intel/composer_xe_2013_sp1.0.080/bin/intel64/ifort
#CC=/opt/intel/composer_xe_2013_sp1.0.080/bin/intel64/icc




###############################################################################

OBJ = cross_kernal.o cuda_utils.o exomol_functions.o GPUManager.o Input.o Timer.o Util.o HITRANStateReader.o ExomolStateReader.o BaseProfile.o VoigtProfile.o DopplerProfile.o StateReader.o MultiGPUManager.o
      # cprio.o

gpu_cross.x:       $(OBJ) main.o
	$(FOR) -o gpu_cross_$(PLAT).x $(OBJ) $(FFLAGS) main.o $(LIB) 

main.o:       main.cu $(OBJ) 
	$(FOR) -c main.cu $(FFLAGS)

cross_kernal.o: cross_kernal.cu cuda_utils.o
	$(FOR) -c cross_kernal.cu $(FFLAGS)

cuda_utils.o: cuda_utils.cu
	$(FOR) -c cuda_utils.cu $(FFLAGS)

HITRANStateReader.o: HITRANStateReader.cpp StateReader.o
	$(FOR) -c HITRANStateReader.cpp $(FFLAGS)

ExomolStateReader.o: ExomolStateReader.cpp StateReader.o
	$(FOR) -c ExomolStateReader.cpp $(FFLAGS)

StateReader.o: StateReader.cpp
	$(FOR) -c StateReader.cpp $(FFLAGS)

BaseProfile.o: BaseProfile.cpp GPUManager.o Input.o MultiGPUManager.o
	$(FOR) -c BaseProfile.cpp $(FFLAGS)

VoigtProfile.o: VoigtProfile.cpp BaseProfile.o Timer.o HITRANStateReader.o ExomolStateReader.o 
	$(FOR) -c VoigtProfile.cpp $(FFLAGS)

DopplerProfile.o: DopplerProfile.cpp BaseProfile.o Timer.o HITRANStateReader.o ExomolStateReader.o
	$(FOR) -c DopplerProfile.cpp $(FFLAGS)

GPUManager.o: GPUManager.cpp
	$(FOR) -c GPUManager.cpp $(FFLAGS)

MultiGPUManager.o: MultiGPUManager.cpp GPUManager.o
	$(FOR) -c MultiGPUManager.cpp $(FFLAGS)
	
Input.o: Input.cpp Util.o
	$(FOR) -c Input.cpp $(FFLAGS)

Util.o: Util.cpp
	$(FOR) -c Util.cpp $(FFLAGS)

Timer.o: Timer.cpp   
	$(FOR) -c Timer.cpp $(FFLAGS)

exomol_functions.o: exomol_functions.cpp Util.o
	$(FOR) -c exomol_functions.cpp $(FFLAGS)

clean:
	rm $(OBJ) *.o main.o

