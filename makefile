goal:   gpu_cross.x



PLAT = GPU




FC = ifort
ICC = icc
ICCFLAGS = -std=c++11 -O3 -xHost
CCFLAGS = -O3 -xHost


ifeq ($(PLAT),GPU)
FOR  = nvcc 
FFLAGS = --ptxas-options=-v -O3 -arch=sm_35 -Xptxas -v -lineinfo --compiler-options "-std=c++0x" -DGPU_ENABLED
CUDA_OBJ = MultiGPUManager.o GPUManager.o cross_kernal.o cuda_utils.o
CUDA_LIB = -lcuda 
MANG_OBJ=MultiGPUManager.o
else
FOR  = $(ICC)
FFLAGS = $(ICCFLAGS)
endif


LIB = -lifcore -limf -lpthread
INC = -I/home/ucapfal/gpu_cross/bz2_compression/

###############################################################################

OBJ =  exomol_functions.o  Input.o Timer.o Util.o HITRANStateReader.o BD_TIPS_2011_v1p0.o ExomolStateReader.o BaseProfile.o VoigtProfile.o BaseManager.o DopplerProfile.o StateReader.o  OpenMPManager.o read_compress_trans.o bzlib.o profiles.o HybridManager.o $(CUDA_OBJ)
      # cprio.o

gpu_cross.x:       $(OBJ) main.o
	$(FOR) -o gpu_cross_$(PLAT).x $(OBJ) main.o /home/ucapfal/bz2/bzip2-1.0.6/libbz2.a $(LIB) $(INC)

main.o:       main.cpp $(OBJ) 
	$(ICC) -c main.cpp $(ICCFLAGS)

cross_kernal.o: cross_kernal.cu cuda_utils.o
	$(FOR) -c cross_kernal.cu --ptxas-options=-v -O3 -arch=sm_35 -Xptxas -v -lineinfo

cuda_utils.o: cuda_utils.cu
	$(FOR) -c cuda_utils.cu --ptxas-options=-v -O3 -arch=sm_35 -Xptxas -v -lineinfo

HITRANStateReader.o: HITRANStateReader.cpp StateReader.o BD_TIPS_2011_v1p0.o 
	$(ICC) -c HITRANStateReader.cpp $(ICCFLAGS)

ExomolStateReader.o: ExomolStateReader.cpp StateReader.o read_compress_trans.o
	$(ICC) -c ExomolStateReader.cpp $(ICCFLAGS)

StateReader.o: StateReader.cpp
	$(FOR) -c StateReader.cpp $(FFLAGS)

BaseProfile.o: BaseProfile.cpp  Input.o  HybridManager.o
	$(FOR) -c BaseProfile.cpp $(FFLAGS)

VoigtProfile.o: VoigtProfile.cpp BaseProfile.o Timer.o HITRANStateReader.o ExomolStateReader.o  
	$(FOR) -c VoigtProfile.cpp $(FFLAGS)

DopplerProfile.o: DopplerProfile.cpp BaseProfile.o Timer.o HITRANStateReader.o ExomolStateReader.o
	$(FOR) -c DopplerProfile.cpp $(FFLAGS)

GPUManager.o: GPUManager.cpp BaseManager.o
	$(FOR) -c GPUManager.cpp $(FFLAGS)

HybridManager.o: HybridManager.cpp BaseManager.o OpenMPManager.o $(MANG_OBJ)
	$(FOR) -c HybridManager.cpp $(FFLAGS)

BaseManager.o: BaseManager.cpp
	$(FOR) -c BaseManager.cpp $(FFLAGS)

BD_TIPS_2011_v1p0.o: 
	ifort -c ./HITRAN_files/BD_TIPS_2011_v1p0.for -O3

MultiGPUManager.o: MultiGPUManager.cpp GPUManager.o BaseManager.o
	$(FOR) -c MultiGPUManager.cpp $(FFLAGS)
	
Input.o: Input.cpp Util.o
	$(FOR) -c Input.cpp $(FFLAGS)

Util.o: Util.cpp
	$(FOR) -c Util.cpp $(FFLAGS)

Timer.o: Timer.cpp   
	$(FOR) -c Timer.cpp $(FFLAGS)

OpenMPManager.o: OpenMPManager.cpp BaseManager.o profiles.o
	$(ICC) -c OpenMPManager.cpp $(ICCFLAGS)

profiles.o: profiles.cpp
	$(ICC) -c profiles.cpp $(CCFLAGS) 

exomol_functions.o: exomol_functions.cpp Util.o
	$(FOR) -c exomol_functions.cpp $(FFLAGS)

bzlib.o: bz2_compression/bzlib.c
	$(ICC) -c bz2_compression/bzlib.c  $(CCFLAGS) 

read_compress_trans.o: bz2_compression/read_compress_trans.c bzlib.o
	$(ICC) -c bz2_compression/read_compress_trans.c $(CCFLAGS) 
clean:
	rm $(OBJ) *.o main.o

