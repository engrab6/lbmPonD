#ifndef _PARAMS_H
#define _PARAMS_H
#include <stdio.h>
//#include <cuda_fp16.h>
#include <cuda_runtime.h>
#include "cuda_math.h"
#include <assert.h>
#include <vector>
#include <string>
#include <typeinfo>
#define STATIC_ASSERT( x ) typedef char __STATIC_ASSERT__[( x )?1:-1]

#define USE_DOUBLE

typedef unsigned short ttype;
//typedef char ttype;
#ifdef USE_FLOAT
#define FTYPESIZE 4
typedef float ftype;
typedef float2 ftype2;
typedef float3 ftype3;
typedef float4 ftype4;
template<typename T1,typename T2>                         __host__ __device__ ftype2 make_ftype2(const T1& f1, const T2& f2) { return make_float2(f1,f2); }
template<typename T1,typename T2,typename T3>             __host__ __device__ ftype3 make_ftype3(const T1& f1, const T2& f2, const T3& f3) { return make_float3(f1,f2,f3); }
template<typename T1,typename T2,typename T3,typename T4> __host__ __device__ ftype4 make_ftype4(const T1& f1, const T2& f2, const T3& f3, const T4& f4) { return make_float4(f1,f2,f3,f4); }
template<typename T> __host__ __device__ ftype3 make_ftype3(const T& f) { return make_float3(f); }
typedef float fptype;
typedef float2 fptype2;
typedef float4 fptype4;
template<typename T1,typename T2> __host__ __device__ fptype2 make_fptype2(const T1& f1, const T2& f2) { return make_float2(f1,f2); }
template<typename T> __host__ __device__ float2 fp22ftype2(const T& f) { return f; }
template<typename T> __host__ __device__ float ftype2fptype(const T& f) { return f; }
template<typename T> __host__ __device__ float fptype2ftype(const T& f) { return f; }
#elif defined USE_DOUBLE
#include "cuda_math_double.h"
#define FTYPESIZE 8
typedef double ftype;
typedef double2 ftype2;
typedef double3 ftype3;
typedef double4 ftype4;
template<typename T1,typename T2>                         __host__ __device__ ftype2 make_ftype2(const T1& f1, const T2& f2) { return make_double2(f1,f2); }
template<typename T1,typename T2,typename T3>             __host__ __device__ ftype3 make_ftype3(const T1& f1, const T2& f2, const T3& f3) { return make_double3(f1,f2,f3); }
template<typename T1,typename T2,typename T3,typename T4> __host__ __device__ ftype4 make_ftype4(const T1& f1, const T2& f2, const T3& f3, const T4& f4) { return make_double4(f1,f2,f3,f4); }
template<typename T> __host__ __device__ ftype3 make_ftype3(const T& f) { return make_double3(f); }
typedef double fptype;
typedef double2 fptype2;
typedef double4 fptype4;
template<typename T1,typename T2> __host__ __device__ fptype2 make_fptype2(const T1& f1, const T2& f2) { return make_double2(f1,f2); }
template<typename T> __host__ __device__ double2 fp22ftype2(const T& f) { return f; }
template<typename T> __host__ __device__ double ftype2fptype(const T& f) { return f; }
template<typename T> __host__ __device__ double fptype2ftype(const T& f) { return f; }
#elif defined USE_HALF
typedef half fptype;
typedef half2 fptype2;
typedef float2 fptype4;
template<typename T1,typename T2> __host__ __device__ fptype2 make_fptype2(const T1& f1, const T2& f2) { return __float22half2_rn(make_float2(f1,f2)); }
template<typename T> __host__ __device__ float2 fp22ftype2(const T& f) { return __half22float2(f); }
template<typename T> __host__ __device__ half ftype2fptype(const T& f) { return __float2half(f); }
template<typename T> __host__ __device__ float fptype2ftype(const T& f) { return __half2float(f); }
#endif

#ifdef USE_UVM
#define GPUALLOC cudaMallocManaged
#define ATTACH ,cudaMemAttachGlobal
#else 
#define GPUALLOC cudaMalloc
#define ATTACH  
#endif

#include "consts.h"

#define DEBUG_PRINT(debug) ;//printf debug;
#define DEBUG_MPI(debug)   ;//printf debug;
#define CHECK_ERROR(err) CheckError( err, __FILE__,__LINE__)
static bool CheckError( cudaError_t err, const char *file, int line) {
  if(err!=cudaSuccess){
    fprintf(stderr, "%s in %s at line %d\n", cudaGetErrorString(err), file, line);
    exit(EXIT_FAILURE);
  }
  return 0;
}
__device__ static double datomicAdd(double* address, double val) {
    unsigned long long int* address_as_ull = (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
      assumed = old;
      old = atomicCAS(address_as_ull, assumed, 
                      __double_as_longlong(val + 
                      __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ static double atomicAdd(double* address, double val) { return datomicAdd(address,val); }
#endif

template<class T> __host__ __device__ inline bool operator==(const T& val1, const T& val2) {
  return val1.x==val2.x && val1.y==val2.y && val1.z==val2.z;
}
template <typename T> inline __host__ __device__ int sign(T val) {
  return (T(0) < val) - (val < T(0));
}
template<class T> __host__ __device__ bool equal(const T& val1, const T& val2) {
  return val1==val2;
}
template<> __host__ __device__ inline bool equal(const int3& val1, const int3& val2) {
  return (val1.x==val2.x && val1.y==val2.y && val1.z==val2.z);
}
__host__ __device__ inline ftype sqrts(const ftype& val) {
  return sqrt(fabs(val))*sign(val);
}
class cuTimer {
  cudaEvent_t tstart,tend;
  cudaStream_t st;
  float diftime;
  public:
  cuTimer(const cudaStream_t& stream=0): diftime(0) {
    CHECK_ERROR( cudaEventCreate(&tstart) ); 
    CHECK_ERROR( cudaEventCreate(&tend  ) );
    CHECK_ERROR( cudaEventRecord(tstart,stream) ); st=stream;
  }
  ~cuTimer(){
    CHECK_ERROR( cudaEventDestroy(tstart) );
    CHECK_ERROR( cudaEventDestroy(tend) );
  }
  float gettime(){
    CHECK_ERROR( cudaEventRecord(tend,st) );
    CHECK_ERROR( cudaEventSynchronize(tend) );
    CHECK_ERROR( cudaEventElapsedTime(&diftime, tstart,tend) ); 
    return diftime;
  }
};

static const int CudaDevs=1;
template<class Ph, class Pd> static void copy2dev(Ph &hostP, Pd &devP) {
  if(CudaDevs>1) {
    int curdev; CHECK_ERROR( cudaGetDevice(&curdev) );
    for(int i=0; i<CudaDevs; i++) {
      CHECK_ERROR( cudaSetDevice(i) );
      CHECK_ERROR( cudaMemcpyToSymbol(devP, &hostP, sizeof(Pd)) );
    }
    CHECK_ERROR( cudaSetDevice(curdev) );
  }
  else CHECK_ERROR( cudaMemcpyToSymbol(devP, &hostP, sizeof(Pd)) );
}
__global__ static void dev_malloc(void* &ptr, size_t size) { ptr=malloc(size); }
__global__ static void dev_free  (void*  ptr) { free(ptr); }
inline void malloc_byGPU(void* &ptr, size_t size) {
  dev_malloc<<<1,1>>>(ptr,size);
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
}
inline void free_byGPU(void* ptr) {
  dev_free<<<1,1>>>(ptr);
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
}
template<class T> void whose_pointer(T* p) {
  printf("Pointer %p: \t",p); fflush(stdout);
  cudaPointerAttributes ptrAt; CHECK_ERROR(cudaPointerGetAttributes(&ptrAt, p));
  if(ptrAt.memoryType==cudaMemoryTypeHost  ) printf("|memoryType: Host     \t");
  if(ptrAt.memoryType==cudaMemoryTypeDevice) printf("|memoryType: Device %d\t", ptrAt.device);
  printf("hostPtr=%p devPtr=%p\n", ptrAt.hostPointer, ptrAt.devicePointer);
}
__device__ inline int ptxAtomicLoadInt(int* global) {
  int local;
  asm("ld.global.cg.s32 %0, [%1];" :"=r"(local) : "l"(global));
//   if(typeid(T)==typeid(         int     )) asm("ld.global.cg.s32 %0, [%1];" :"=r"(local) : "l"(global));
//   if(typeid(T)==typeid(unsigned int     )) asm("ld.global.cg.u32 %0, [%1];" :"=r"(local) : "l"(global));
//   if(typeid(T)==typeid(         long int)) asm("ld.global.cg.s64 %0, [%1];" :"=l"(local) : "l"(global));
//   if(typeid(T)==typeid(unsigned long int)) asm("ld.global.cg.u64 %0, [%1];" :"=l"(local) : "l"(global));
//   if(typeid(T)==typeid(float            )) asm("ld.global.cg.f32 %0, [%1];" :"=f"(local) : "l"(global));
//   if(typeid(T)==typeid(double           )) asm("ld.global.cg.f64 %0, [%1];" :"=d"(local) : "l"(global));
  return local;
}

__host__ __device__ inline static long morton_zip(unsigned long int x, unsigned long int y, unsigned long int z) {
  //------ code based on https://graphics.stanford.edu/~seander/bithacks.html
  auto drawl = [] (unsigned long int& v) { 
    static const unsigned long int B[] = {0x249249249249UL, 0xC30C30C30C3UL, 0xF00F00F00FUL, 0xFF0000FFUL};
    static const unsigned int S[] = {2, 4, 8, 16};
    v = (v | (v << S[3])) & B[3];
    v = (v | (v << S[2])) & B[2];
    v = (v | (v << S[1])) & B[1];
    v = (v | (v << S[0])) & B[0];
  };
   // Interleave lower 16 bits of x and y, so the bits of x
   // are in the even positions and bits from y in the odd;
   // z gets the resulting 32-bit Morton Number.  
   // x and y must initially be less than 65536.
  drawl(x);
  drawl(y);
  drawl(z);
  return x | (y << 1) | (z << 2);
}

#include "im2D.h"
#include "im3D.hpp"
extern image2D im2D;
extern float calcTime, calcPerf; extern int TimeStep;
extern bool recalc_at_once;
extern const char* FuncStr[];

#endif //_PARAMS_H
