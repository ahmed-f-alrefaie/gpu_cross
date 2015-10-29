#ifndef CROSS_KERNAL_H
#define CROSS_KERNAL_H

#include "cross_structs.cuh"
#include <cuda_runtime_api.h>
#include <cuda.h>
#include <cmath>

#ifdef KEPLER
#define BLOCK_SIZE 1024
#else
#define BLOCK_SIZE 576
#endif
#define SHARED_SIZE 256
#define VOIGT_BLOCK 128
#define VOIGT_SHARED_SIZE 128

__host__ void copy_intensity_info(cross_section_data* cross_inf);

//__global__ void device_compute_cross_section(const double* g_freq, double* g_cs,const double* g_energies,const int* g_gns,const int* g_jF,const double* g_nu,const double* g_aif, const int N,const int N_ener,const int start_idx);

__global__ void device_compute_cross_section(const double*  __restrict__ g_freq, double* g_cs,const double*  __restrict__ g_energies,const int*  __restrict__ g_gns,const int*  __restrict__ g_jF,const double*  __restrict__ g_nu,const double*  __restrict__ g_aif, const int N,const int N_ener,const int start_idx);

__global__ void device_compute_cross_section_abscoef(const double*  __restrict__ g_freq, double* g_cs,const double*  __restrict__ g_energies,const int*  __restrict__ g_gns,const double*  __restrict__ g_nu,const double*  __restrict__ g_aif, const int N,const int N_ener,const int start_idx);

__global__ void device_compute_cross_section_warp_abscoef(const double*  __restrict__ g_freq, double* g_cs,const double*  __restrict__ g_energies,const int*  __restrict__ g_gns,const double*  __restrict__ g_nu,const double*  __restrict__ g_aif, const int N,const int N_ener,const int start_idx);

__global__ void device_compute_cross_section_stepone(double* g_energies,const int*  __restrict__ g_gns,const double*  __restrict__ g_nu,const double*  __restrict__ g_aif, const int N_ener);
__global__ void device_compute_cross_section_steptwo(const double*  __restrict__ g_freq, double* g_cs,const double*  __restrict__ g_nu,const double*  __restrict__ g_abscoef,const int N,const int N_ener,const int start_idx);

__host__ void execute_two_step_kernal(double* g_freq, double* g_intens, double* g_energies, double* g_nu, int* g_gns,double* g_aif, int Npoints,int N_ener);

__global__ void device_compute_cross_section_steptwo_block(const double* g_freq, double* g_cs,const double* g_nu,const double* g_abscoef,const int N,const int N_ener,const int start_idx);

__host__ void execute_two_step_kernal_block(double* g_freq, double* g_intens, double* g_energies, double* g_nu, int* g_gns,double* g_aif, int Npoints,int N_ener,int start_idx);

//VOIGT

__global__ void device_compute_cross_section_voigt_stepone(double* g_energies,const int*  __restrict__ g_gns,const double*  __restrict__ g_nu,const double*  __restrict__ g_aif, const int N_ener);
__global__ void device_compute_cross_section_voigt_stepone(double* g_energies,const int*  g_gns,const double*  g_nu,const double*  g_aif,double*  g_gamma,double*  g_n, const int N_ener);
__global__ void device_compute_cross_section_voigt_steptwo_block(const double*  g_freq, double* g_cs,const double*   g_nu,const double*  g_abscoef,const double*  g_gamma,const int N,const int N_ener,const int start_idx);
__global__ void device_compute_cross_section_voigt_steptwo_block(const double*  g_freq, double* g_cs,const double*   g_nu,const double*  g_abscoef,const int N,const int N_ener,const int start_idx);

__host__ void execute_two_step_kernal_voigt_block(double* g_freq, double* g_intens, double* g_energies, double* g_nu, int* g_gns,double* g_aif, int Npoints,int N_ener,int start_idx);
__host__ void execute_two_step_kernal_voigt_block(double* g_freq, double* g_intens, double* g_energies, double* g_nu, int* g_gns,double* g_aif,double* g_gamma,double* g_n ,int Npoints,int N_ener,int start_idx);
__host__ void execute_two_step_kernal_voigt(double* g_freq, double* g_intens, double* g_energies, double* g_nu, int* g_gns,double* g_aif, int Npoints,int N_ener,const int start_idx);

__global__ void device_compute_cross_section_voigt_steptwo(const double* g_freq, double* g_cs,const double*  g_nu,const double*  g_abscoef,const int N,const int N_ener,const int start_idx);


//////DOPPLER////////////
__host__ void execute_two_step_kernal_doppler_block(double* g_freq, double* g_intens, double* g_energies, double* g_nu, int* g_gns,double* g_aif, int Npoints,int N_ener,int start_idx);
__global__ void device_compute_cross_section_doppler_stepone(double* g_energies,const int*  __restrict__ g_gns,const double*  __restrict__ g_nu,const double*  __restrict__ g_aif, const int N_ener);
__global__ void device_compute_cross_section_doppler_steptwo_block(const double*  g_freq, double* g_cs,const double*   g_nu,const double*  g_abscoef,const int N,const int N_ener,const int start_idx);



#endif
