#pragma once
#include <cuda.h>
#include <curand_kernel.h>
#include "params.h"

const int WorldLine=1;
const long int Gblock=1;
const int Gsize=1;
const long int Rank=4;
const long int LSizeX=32;//32;
const long int LSizeY=1;//32;//4;
const long int LSizeZ=1;//32;//16*4;//4;
const long int Nt=1;//1<<Rank;
// const long int tfNx=((1<<Rank)+Nt-1)/Gblock;
// const long int tfNy=((1<<Rank)+Nt-1)/Gblock;
// const long int tfNz=((1<<Rank)+Nt-1)/Gblock;
static constexpr uint3 brick=(const uint3){4u,1u,1u};
static constexpr int MultiplyDeBruijnBitPosition[32] = { 0,9,1,10,13,21,2,29,11,14,16,18,22,25,3,30,8,12,20,28,15,17,24,7,19,27,23,6,26,5,4,31 };
constexpr __host__ __device__ int shiftOR(int v, int shift) { return v|v>>shift; }
constexpr __host__ __device__ int MSB(int v) {
  return MultiplyDeBruijnBitPosition[(uint32_t)(shiftOR( shiftOR( shiftOR( shiftOR( shiftOR(v,1), 2 ), 4 ), 8), 16 ) * 0x07C4ACDDU) >> 27];
}
// const long int tfNx = 1<<(MSB( ((1<<Rank)+Nt-1)/Gblock-1 )+1);
// const long int tfNy = 1<<(MSB( ((1<<Rank)+Nt-1)/Gblock-1 )+1);
// const long int tfNz = 1<<(MSB( ((1<<Rank)+Nt-1)/Gblock-1 )+1);
const long int tfNx = 1<<(MSB( ((1<<Rank)+(WorldLine*Nt-1)/brick.x+1)/Gblock-1 )+1);
const long int tfNy = 1<<(MSB( ((1<<Rank)+(WorldLine*Nt-1)/brick.y+1)/Gblock-1 )+1);
const long int tfNz = 1<<(MSB( ((1<<Rank)+(WorldLine*Nt-1)/brick.z+1)/Gblock-1 )+1);

struct CellP{
  ftype4 f[5];
  ftype4 uT;
};
struct Group;
struct Cell{
  ftype f[Qn];
  #ifdef DDF
  ftype h[Qn];
  #endif
  #ifdef NRMESH
  ftype3 p;
  #endif
  ftype rho;
  ftype4 uT;
  static long __host__ __device__ ind_zip(long int x, long int y, long int z);
  void __host__ __device__ set(const ftype _r, const ftype4 _u);
  void __host__ __device__ set_base(const ftype _r, const ftype4 _u);
  void __device__ calc_eq(ftype feq[Qn]);
  void __device__ collision();
  void __device__ fast_collision();
  void __device__ gather(const long index) {}
  void __device__ gather(const int x, const int y, const int z);
  void __device__ scatter(const int x, const int y, const int z);
  void __device__ gather_group(const int x, const int y, const int z, Group* g1=0, Group* g2=0);
  void __device__ scatter_group(const int x, const int y, const int z);
  void __device__ load   (const long int x, const long int y, const long int z);
  void __device__ load_st(const long int x, const long int y, const long int z);
  void __device__ save   (const long int x, const long int y, const long int z);
  void __device__ save_ld(const long int x, const long int y, const long int z);
  void __device__ save_12(const long int x, const long int y, const long int z);
  ftype2 __device__ transfer_gauge(ftype4 base, int inew, int ix=0, ftype M4val=0);
  void __host__ __device__ calcMoments( ftype3 vi[Qn], ftype ux=0 );
  void __device__ mirrX();
  void __device__ mirrY();
  void __device__ mirrZ();
};
struct CellsSoA{
  ftype* f[Qn];
  ftype* rho;
  static long __host__ __device__ ind_zip(long int x, long int y, long int z);
  void __host__ __device__ set(const long int index, const ftype _r, const ftype4 _u);
  void __device__ load (Cell& c, const long int index);
  void __device__ store(Cell& c, const long int index);
  void __device__ gather (Cell& c, const int x, const int y, const int z);
  void __device__ scatter(Cell& c, const int x, const int y, const int z);
};
struct Group{
  ftype f[Qn][Gsize*Gsize*Gsize];
  ftype h[Qn][Gsize*Gsize*Gsize];
  ftype3 p[Gsize*Gsize*Gsize];
  ftype rho[Gsize*Gsize*Gsize];
  ftype4 uT[Gsize*Gsize*Gsize];
  static long __host__ __device__ ind_zip(long int x, long int y, long int z) { if(Gsize==1) return Cell::ind_zip(x,y,z); else return morton_zip(x,y,z); };
  static int2 __host__ __device__ ind_conv(long int x, long int y, long int z);
  void __host__ __device__ pack  (Cell& c, int cindex);
  void __host__ __device__ unpack(Cell& c, int cindex);
};

struct Cinfo{
  ftype ctype;
  int niter;
  int set,setnew;
  static long __host__ __device__ ind_zip(long int x, long int y, long int z) { return Cell::ind_zip(x,y,z); };
};

struct PhysPars {
  ftype dtau;
  void set() { dtau=1.0/0.9;/*0.51*/; }
};
struct LBMParams{
  long int Nx,Ny,Nz;
  unsigned int nFunc;
  int iStep;
  double curtime;
  double dt;
  int* NiterMax;
  ftype* mass, *enrg;
  ftype3* moment;
  cudaArray* data;
  cudaArray* fdata;
  cudaArray* gdata;
  cudaSurfaceObject_t texSdata;
  cudaSurfaceObject_t texFdata;
  cudaSurfaceObject_t texGdata;
  Cell* cells;
  Cinfo* cinfo;
  Group* groups;
  Group* loadGroups;
  Group* storeGroups;
  CellsSoA csoa;
  PhysPars pc;
};

struct LBMParamsHost: public LBMParams {
  Arr3D_pars arr4im, arr4surf;
  unsigned int MaxFunc;
  void set();
  void reset();
  void reset_im() {
    int ndevs=0; int curdev;
    CHECK_ERROR( cudaGetDevice(&curdev) );
    CHECK_ERROR( cudaGetDeviceCount(&ndevs) );
    //if(ndevs>1) CHECK_ERROR(cudaSetDevice(1));
    arr4im.reset(Nx,Ny,Nz);
    arr4im.BufSize = sizeof(float)*Nx*Ny*Nz;
    CHECK_ERROR( cudaMallocHost((void**) (&arr4im.Arr3Dbuf), arr4im.BufSize) );
    CHECK_ERROR( cudaMemset(arr4im.Arr3Dbuf, 0, arr4im.BufSize) );
    arr4im.inGPUmem = true;
    arr4surf.inGPUmem = true;
    CHECK_ERROR(cudaSetDevice(curdev));
  }
  void clear() {
    //cudaFree (g1); cudaFree (g2);
    cudaFree (arr4im.Arr3Dbuf);
    arr4im.clear();
  }
  void checkSizes() { }
  void swapLS() { Group* tmp=storeGroups; storeGroups=loadGroups; loadGroups=tmp; }
};
extern LBMParamsHost parsHost;
extern __constant__ LBMParams pars;

inline long Cell::ind_zip(long int x, long int y, long int z) {
//  return morton_zip(x,y,z);
//   long int xb=x/brick.x, yb=y/brick.y, zb=z/brick.z;
//   #ifdef __CUDA_ARCH__
//   return (xb*brick.x*brick.y*brick.z+yb*pars.Ny    *brick.x*brick.z+zb*pars.Nx*pars.Ny*brick.z)+x%brick.x+y%brick.y*brick.x+z%brick.z*brick.x*brick.y;
//   #else
//   return (xb*brick.x*brick.y*brick.z+yb*parsHost.Ny*brick.x*brick.z+zb*parsHost.Nx*parsHost.Ny*brick.z)+x%brick.x+y%brick.y*brick.x+z%brick.z*brick.x*brick.y;
//   #endif
  #ifdef __CUDA_ARCH__
  return x+y*pars.Nx+z*pars.Nx*pars.Ny;
  #else
  return x+y*parsHost.Nx+z*parsHost.Nx*parsHost.Ny;
  #endif
}
inline long CellsSoA::ind_zip(long int x, long int y, long int z) {
  return Cell::ind_zip(x,y,z);
}
inline void __device__ CellsSoA::load (Cell& c, const long int index){
  for(int i=0; i<Qn; i++) c.f[i]=f[i][index];
  c.rho=rho[index];
}
inline void __device__ CellsSoA::store(Cell& c, const long int index){
  for(int i=0; i<Qn; i++) f[i][index]=c.f[i];
  rho[index]=c.rho;
}
inline void __device__ CellsSoA::gather(Cell& c, const int cx, const int cy, const int cz){
  for(int i=Qn-1; i>=0; i--) {
    int3 crd = make_int3(cx,cy,cz);
    int3 ei=e[reverseXYZ[i]];
    int inear=reverseXYZ[i];
    crd+= e[i];
    if(crd.x<0 || crd.y<0 || crd.z<0 || crd.x>=pars.Nx || crd.y>=pars.Ny || crd.z>=pars.Nz) {
      crd = make_int3(cx,cy,cz);
      inear=i;
    }
    const long index = ind_zip(crd.x, crd.y, crd.z);
    c.f[i] = f[inear][index];
  }
}
inline void __device__ CellsSoA::scatter(Cell& c, const int cx, const int cy, const int cz){
  for(int i=Qn-1; i>=0; i--) {
    int3 crd = make_int3(cx,cy,cz);
    int3 ei=e[reverseXYZ[i]];
    int inear=reverseXYZ[i];
    crd+= e[i];
    if(crd.x<0 || crd.y<0 || crd.z<0 || crd.x>=pars.Nx || crd.y>=pars.Ny || crd.z>=pars.Nz) {
      crd = make_int3(cx,cy,cz);
      inear=i;
    }
    const long index = ind_zip(crd.x, crd.y, crd.z);
    f[inear][index] = c.f[i];
  }
}
inline void __device__ Cell::gather(const int cx, const int cy, const int cz){
  for(int i=Qn-1; i>=0; i--) {
    int3 crd = make_int3(cx,cy,cz);
    int3 ei=e[reverseXYZ[i]];
    int inear=reverseXYZ[i];
    crd+= e[i];
    if(crd.x<0 || crd.y<0 || crd.z<0 || crd.x>=pars.Nx || crd.y>=pars.Ny || crd.z>=pars.Nz) {
      crd = make_int3(cx,cy,cz);
      inear=i;
    }
    const long index = ind_zip(crd.x, crd.y, crd.z);
    f[i] = pars.cells[index].f[inear];
  }
}
inline void __device__ Cell::scatter(const int cx, const int cy, const int cz){
  for(int i=Qn-1; i>=0; i--) {
    int3 crd = make_int3(cx,cy,cz);
    int3 ei=e[reverseXYZ[i]];
    int inear=reverseXYZ[i];
    crd+= e[i];
    if(crd.x<0 || crd.y<0 || crd.z<0 || crd.x>=pars.Nx || crd.y>=pars.Ny || crd.z>=pars.Nz) {
      crd = make_int3(cx,cy,cz);
      inear=i;
    }
    const long index = ind_zip(crd.x, crd.y, crd.z);
    pars.cells[index].f[inear] = f[i];
  }
}
inline void __device__ Cell::gather_group(const int cx, const int cy, const int cz, Group* gptr1, Group* gptr2){
 for(int i=0; i<Qn; i++) {
    int3 crd = make_int3(cx,cy,cz);
    int3 ei=e[reverseXYZ[i]];
    int inear=reverseXYZ[i];
    crd+= e[i];
    if(crd.x<0 || crd.y<0 || crd.z<0 || crd.x>=pars.Nx || crd.y>=pars.Ny || crd.z>=pars.Nz) {
      crd = make_int3(cx,cy,cz);
      inear=i;
    }
    const int2 gcind = Group::ind_conv(crd.x, crd.y, crd.z);
    const int gindex=gcind.x, cindex=gcind.y;
    f[i] = pars.loadGroups[gindex].f[inear][cindex];
  }
}
inline void __device__ Cell::scatter_group(const int cx, const int cy, const int cz){
  for(int i=0; i<Qn; i++) {
    int3 crd = make_int3(cx,cy,cz);
    int3 ei=e[reverseXYZ[i]];
    register int inear=reverseXYZ[i];
    crd+= e[i];
    if(crd.x<0 || crd.y<0 || crd.z<0 || crd.x>=pars.Nx || crd.y>=pars.Ny || crd.z>=pars.Nz) {
      crd = make_int3(cx,cy,cz);
      inear=i;
    }
    const int2 gcind = Group::ind_conv(crd.x, crd.y, crd.z);
    const int gindex=gcind.x, cindex=gcind.y;
    pars.storeGroups[gindex].f[inear][cindex] = f[i];
  }
}
inline int2 Group::ind_conv(long int x, long int y, long int z){
  const int gindex = ind_zip(x/Gsize, y/Gsize, z/Gsize);
  const int cindex = x%Gsize+y%Gsize*Gsize+z%Gsize*Gsize*Gsize;
  return make_int2(gindex,cindex);
}
inline void __host__ __device__ Group::pack  (Cell& c, int cindex){
  for(int i=0; i<Qn; i++) f[i][cindex]=c.f[i];
  #ifdef DDF
  for(int i=0; i<Qn; i++) h[i][cindex]=c.h[i];
  #endif
  #ifdef NRMESH
  p[cindex]=c.p;
  #endif
  rho[cindex]=c.rho; uT[cindex]=c.uT;
}
inline void __host__ __device__ Group::unpack(Cell& c, int cindex){
  for(int i=0; i<Qn; i++) c.f[i]=f[i][cindex];
  #ifdef DDF
  for(int i=0; i<Qn; i++) c.h[i]=h[i][cindex];
  #endif
  #ifdef NRMESH
  c.p=p[cindex];
  #endif
  c.rho=rho[cindex]; c.uT=uT[cindex];
}
inline void __device__ Cell::load(const long int x, const long int y, const long int z){
  const int2 gcind = Group::ind_conv(x, y, z);
  const int gindex=gcind.x, cindex=gcind.y;
  pars.loadGroups[gindex].unpack(*this, cindex);
}
inline void __device__ Cell::load_st(const long int x, const long int y, const long int z){
  const int2 gcind = Group::ind_conv(x, y, z);
  const int gindex=gcind.x, cindex=gcind.y;
  pars.storeGroups[gindex].unpack(*this, cindex);
}
inline void __device__ Cell::save(const long int x, const long int y, const long int z){
  const int2 gcind = Group::ind_conv(x, y, z);
  const int gindex=gcind.x, cindex=gcind.y;
  pars.storeGroups[gindex].pack(*this, cindex);
}
inline void __device__ Cell::save_ld(const long int x, const long int y, const long int z){
  const int2 gcind = Group::ind_conv(x, y, z);
  const int gindex=gcind.x, cindex=gcind.y;
  pars.loadGroups[gindex].pack(*this, cindex);
}
inline void __device__ Cell::save_12(const long int x, const long int y, const long int z){
  const int2 gcind = Group::ind_conv(x, y, z);
  const int gindex=gcind.x, cindex=gcind.y;
  pars.loadGroups [gindex].pack(*this, cindex);
  pars.storeGroups[gindex].pack(*this, cindex);
}
inline void __host__ __device__ Cell::calcMoments(ftype3 vi[Qn], ftype ux){
  rho=0; uT= make_ftype4(0,0,0,0);
  ftype M3=0,Ma=0,Mb=0,Mc=0;
  /*Ma= f[0]*vi[0].x+f[1]*vi[1].x+f[2]*vi[2].x+f[3]*(ux)+f[4]*(ux);
  Mb= f[0]*vi[0].x+f[3]*vi[3].x+f[4]*vi[4].x+f[1]*(ux)+f[2]*(ux);
  ftype corrcoeff=1;
  if(f[1]-f[2]<0 && f[3]-f[4]<0) corrcoeff=-1;
  Mc= f[0]*vi[0].x+f[1]*ux+f[2]*ux+f[3]*ux+f[4]*ux+corrcoeff*(fabs(f[1]*(vi[1].x-ux)+f[2]*(vi[2].x-ux))+fabs(f[3]*(vi[3].x-ux)+f[4]*(vi[4].x-ux)));*/
  for(int i=0; i<Qn; i++) {
    rho+=f[i];
    uT+= make_ftype4(f[i]*vi[i].x, f[i]*vi[i].y, f[i]*vi[i].z, f[i]*dot(vi[i],vi[i]));
    M3+= f[i]*vi[i].x*vi[i].x*vi[i].x; 
  }
  //uT.x=Mc;
  if(rho!=0) uT/= rho;
  uT.w = (uT.w-uT.x*uT.x-uT.y*uT.y-uT.z*uT.z)*(1./Dim);
  ftype ufromM3 = M3/(rho*(uT.x*uT.x+3*uT.w));
  ftype TfromM3 = (M3/(rho*uT.x)-uT.x*uT.x)/3;
  //if(fabs(uT.x)<1e-7) TfromM3=uT.w;
  //uT.w=(uT.w+TfromM3)*0.5;
  //uT.x = (uT.x+ufromM3)*0.5;
#ifdef DDF
  uT.w= 0;
  for(int i=0; i<Qn; i++) uT.w+=h[i];
  if(rho!=0) uT.w/= rho;
  uT.w = (uT.w-uT.x*uT.x-uT.y*uT.y-uT.z*uT.z)*(1./Dim);
#endif
}
inline ftype2 __device__ Cell::transfer_gauge(ftype4 base, int inew, int ix, ftype M4){
  // from 0 -> to 1
  const ftype Tbase=TLat; //1./3;
  ftype3 u0 = make_ftype3(uT.x,uT.y,uT.z);
  ftype T0=uT.w;
  ftype3 u1 = make_ftype3(base.x,base.y,base.z);
  ftype T1=base.w;
  ftype fnew = 0,hnew=0;
  ftype wk = 1./(2*T1),wl=wk,wm=wk;
  #ifdef D1Q5
  ftype e1_fr=ef[1].x*sqrts(T0/Tbase);
  ftype e2_fr=ef[3].x*sqrts(T0/Tbase);
  ftype e1_to=ef[1].x*sqrts(T1/Tbase);
  ftype e2_to=ef[3].x*sqrts(T1/Tbase);
  ftype e_fr[Qn], e_to[Qn];
  for(int i=0; i<Qn; i++) {
    e_fr[i]=ef[i].x*sqrts(T0/Tbase);
    e_to[i]=ef[i].x*sqrts(T1/Tbase);
  }
  if(ef[inew].x==0       ) wk = 1./(e1_to*e1_to*e2_to*e2_to); else
  if(abs(ef[inew].x)==ec1) wk = 1./(2*e1_to*e1_to*(e1_to*e1_to-e2_to*e2_to));
  else                     wk = 1./(2*e2_to*e2_to*(e2_to*e2_to-e1_to*e1_to));
  ftype M4moment=0;
  //if(ix==1) M4moment=3*T1*T1+6*T1*u1.x*u1.x + u1.x*u1.x*u1.x*u1.x;
  //else M4moment=3*T0*T0+6*T0*u0.x*u0.x + u0.x*u0.x*u0.x*u0.x;
  fnew = wk*M4moment;
  for(int iold=0; iold<Qn; iold++) {
    ftype K0,K1,K2,K3,K4;
    K4=1;
    K3 = -4*e_fr[iold]-e_to[inew];
    if(inew==0) {
      K2 = 6*e_fr[iold]*e_fr[iold]-(e1_to*e1_to+e2_to*e2_to);
      K1 = 2*e_fr[iold]*(e1_to*e1_to+e2_to*e2_to-2*e_fr[iold]*e_fr[iold]);
      K0 = (e1_to*e1_to-e_fr[iold]*e_fr[iold])*(e2_to*e2_to-e_fr[iold]*e_fr[iold]);
    } else {
      K2 = 6*e_fr[iold]*e_fr[iold]+3*e_fr[iold]*e_to[inew]-(e1_to*e1_to+e2_to*e2_to-e_to[inew]*e_to[inew]);
      K1 = (2*e_fr[iold]+e_to[inew])*(e1_to*e1_to+e2_to*e2_to-e_to[inew]*e_to[inew])-e_fr[iold]*e_fr[iold]*(4*e_fr[iold]+3*e_to[inew]);
      K0 = e_fr[iold]*(e_fr[iold]*e_fr[iold]-(e1_to*e1_to+e2_to*e2_to)+e_to[inew]*e_to[inew])*(e_fr[iold]+e_to[inew]);
    }
    const ftype difu=u1.x-u0.x;
    ftype difM4moment = u0.x+e_fr[iold]; difM4moment = difM4moment*difM4moment*difM4moment*difM4moment;
    difM4moment*=0;
    //if(fabs(M4)>0.001) difM4moment*=-0.5;
    //else difM4moment*=0;
    const ftype gx = K4*difu*difu*difu*difu + K3*difu*difu*difu + K2*difu*difu + K1*difu + K0 - difM4moment;
    ftype fadd=0,hadd=0;
    if(Dim==1) fadd=wk*gx*f[iold];
    if(!isnan(fadd) && !isinf(fadd)) fnew+= fadd;
    #ifdef DDF
    if(Dim==1) hadd=wk*gx*h[iold];
    if(!isnan(hadd) && !isinf(hadd)) hnew+= hadd;
    #endif
  }
  #elif defined D3Q125
  const ftype S0 = sqrts(T0/Tbase);
  const ftype S1 = sqrts(T1/Tbase);
  ftype e1_fr=ec1*S0;
  ftype e2_fr=ec2*S0;
  ftype e1_to=ec1*S1;
  ftype e2_to=ec2*S1;
  ftype3 e_fr[Qn], e_to[Qn];
  for(int i=0; i<Qn; i++) {
    e_fr[i]=ef[i]*S0;
    e_to[i]=ef[i]*S1;
  }
  if(ef[inew].x==0       ) wk = 1./(e1_to*e1_to*e2_to*e2_to); else
  if(abs(ef[inew].x)==ec1) wk = 1./(2*e1_to*e1_to*(e1_to*e1_to-e2_to*e2_to));
  else                     wk = 1./(2*e2_to*e2_to*(e2_to*e2_to-e1_to*e1_to));
  if(ef[inew].y==0       ) wl = 1./(e1_to*e1_to*e2_to*e2_to); else
  if(abs(ef[inew].y)==ec1) wl = 1./(2*e1_to*e1_to*(e1_to*e1_to-e2_to*e2_to));
  else                     wl = 1./(2*e2_to*e2_to*(e2_to*e2_to-e1_to*e1_to));
  if(ef[inew].z==0       ) wm = 1./(e1_to*e1_to*e2_to*e2_to); else
  if(abs(ef[inew].z)==ec1) wm = 1./(2*e1_to*e1_to*(e1_to*e1_to-e2_to*e2_to));
  else                     wm = 1./(2*e2_to*e2_to*(e2_to*e2_to-e1_to*e1_to));
  const ftype3 difu=u1-u0;
  const ftype3 difu2=difu*difu;
  const ftype3 difu3=difu*difu*difu;
  const ftype3 difu4=difu*difu*difu*difu;
  const ftype e12to2 = e1_to*e1_to+e2_to*e2_to;
  const ftype3 du=u0-u1;
  int Xindex=0,Yindex=0,Zindex=0;
  if(ef[inew].x==+ec1) Xindex=1; else if(ef[inew].x==-ec1) Xindex=2; else if(ef[inew].x==+ec2) Xindex=3; else if(ef[inew].x==-ec2) Xindex=4;
  if(ef[inew].y==+ec1) Yindex=1; else if(ef[inew].y==-ec1) Yindex=2; else if(ef[inew].y==+ec2) Yindex=3; else if(ef[inew].y==-ec2) Yindex=4;
  if(ef[inew].z==+ec1) Zindex=1; else if(ef[inew].z==-ec1) Zindex=2; else if(ef[inew].z==+ec2) Zindex=3; else if(ef[inew].z==-ec2) Zindex=4;
  for(int iold=0; iold<Qn; iold++) {
    const ftype factX[5] = { (du.x+e_fr[iold].x), (du.x+e_fr[iold].x-e1_to), (du.x+e_fr[iold].x+e1_to), (du.x+e_fr[iold].x-e2_to), (du.x+e_fr[iold].x+e2_to) };
    const ftype factY[5] = { (du.y+e_fr[iold].y), (du.y+e_fr[iold].y-e1_to), (du.y+e_fr[iold].y+e1_to), (du.y+e_fr[iold].y-e2_to), (du.y+e_fr[iold].y+e2_to) };
    const ftype factZ[5] = { (du.z+e_fr[iold].z), (du.z+e_fr[iold].z-e1_to), (du.z+e_fr[iold].z+e1_to), (du.z+e_fr[iold].z-e2_to), (du.z+e_fr[iold].z+e2_to) };

    ftype gx=factX[(Xindex+1)%5]*factX[(Xindex+2)%5]*factX[(Xindex+3)%5]*factX[(Xindex+4)%5];
    ftype gy=factY[(Yindex+1)%5]*factY[(Yindex+2)%5]*factY[(Yindex+3)%5]*factY[(Yindex+4)%5];
    ftype gz=factZ[(Zindex+1)%5]*factZ[(Zindex+2)%5]*factZ[(Zindex+3)%5]*factZ[(Zindex+4)%5];

    /*ftype3 K0,K1,K2,K3,K4;
    ftype3 K0_0,K1_0,K2_0;
    ftype3 K0_1,K1_1,K2_1;
    K4 = make_ftype3(1,1,1);
    K3 = -4*e_fr[iold]-e_to[inew];*/

    /*K2_0 = 6*e_fr[iold]*e_fr[iold]-(e1_to*e1_to+e2_to*e2_to);
    K1_0 = 2*e_fr[iold]*(e1_to*e1_to+e2_to*e2_to-2*e_fr[iold]*e_fr[iold]);
    K0_0 = (e1_to*e1_to-e_fr[iold]*e_fr[iold])*(e2_to*e2_to-e_fr[iold]*e_fr[iold]);
      
    K2_1 = 6*e_fr[iold]*e_fr[iold]+3*e_fr[iold]*e_to[inew]-(e1_to*e1_to+e2_to*e2_to-e_to[inew]*e_to[inew]);
    K1_1 = (2*e_fr[iold]+e_to[inew])*(e1_to*e1_to+e2_to*e2_to-e_to[inew]*e_to[inew])-e_fr[iold]*e_fr[iold]*(4*e_fr[iold]+3*e_to[inew]);
    K0_1 = e_fr[iold]*(e_fr[iold]*e_fr[iold]-(e1_to*e1_to+e2_to*e2_to)+e_to[inew]*e_to[inew])*(e_fr[iold]+e_to[inew]);
                
    if(ef[inew].x==0) { K0.x=K0_0.x; K1.x=K1_0.x; K2.x=K2_0.x; } else { K0.x=K0_1.x; K1.x=K1_1.x; K2.x=K2_1.x; }
    if(ef[inew].y==0) { K0.y=K0_0.y; K1.y=K1_0.y; K2.y=K2_0.y; } else { K0.y=K0_1.y; K1.y=K1_1.y; K2.y=K2_1.y; }
    if(ef[inew].z==0) { K0.z=K0_0.z; K1.z=K1_0.z; K2.z=K2_0.z; } else { K0.z=K0_1.z; K1.z=K1_1.z; K2.z=K2_1.z; }*/
    
//     ftype3 gcoff = (du+e_fr[iold])*(du+e_fr[iold]-e1_to)*(du+e_fr[iold]+e1_to)*(du+e_fr[iold]-e2_to)*(du+e_fr[iold]+e2_to)/(du+e_fr[iold]-e_to[inew]);
//     ftype gxyz=gcoff.x*gcoff.y*gcoff.z;
    /*
    ftype gx = K4.x*difu4.x + K3.x*difu3.x; 
    ftype gy = K4.y*difu4.y + K3.y*difu3.y; 
    ftype gz = K4.z*difu4.z + K3.z*difu3.z;
    if(ef[inew].x==0) {
      ftype efrom = e_fr[iold].x;
      ftype efrom2 = efrom*efrom;
      ftype K0 = (e1_to*e1_to-efrom2)*(e2_to*e2_to-efrom2);
      ftype K1 = 2*efrom*(e12to2-2*efrom2);
      ftype K2 = 6*efrom2-e12to2;
      gx+= K2*difu2.x + K1*difu.x + K0;
    } else {
      ftype efrom = e_fr[iold].x;
      ftype efrom2 = efrom*efrom;
      ftype eto = e_to[inew].x;
      ftype K0 =  efrom*(efrom2-e12to2+eto*eto)*(efrom+eto);
      ftype K1 =  (2*efrom+eto)*(e12to2-eto*eto)-efrom2*(4*efrom+3*eto);
      ftype K2 =  6*efrom2+3*efrom*eto-(e12to2-eto*eto);
      gx+= K2*difu2.x + K1*difu.x + K0;
    }
    if(ef[inew].y==0) {
      ftype efrom = e_fr[iold].y;
      ftype efrom2 = efrom*efrom;
      ftype K0 = (e1_to*e1_to-efrom2)*(e2_to*e2_to-efrom2);
      ftype K1 = 2*efrom*(e12to2-2*efrom2);
      ftype K2 = 6*efrom2-e12to2;
      gy+= K2*difu2.y + K1*difu.y + K0;
    } else {
      ftype efrom = e_fr[iold].y;
      ftype efrom2 = efrom*efrom;
      ftype eto = e_to[inew].y;
      ftype K0 =  efrom*(efrom2-e12to2+eto*eto)*(efrom+eto);
      ftype K1 =  (2*efrom+eto)*(e12to2-eto*eto)-efrom2*(4*efrom+3*eto);
      ftype K2 =  6*efrom2+3*efrom*eto-(e12to2-eto*eto);
      gy+= K2*difu2.y + K1*difu.y + K0;
    }
    if(ef[inew].z==0) {
      ftype efrom = e_fr[iold].z;
      ftype efrom2 = efrom*efrom;
      ftype K0 = (e1_to*e1_to-efrom2)*(e2_to*e2_to-efrom2);
      ftype K1 = 2*efrom*(e12to2-2*efrom2);
      ftype K2 = 6*efrom2-e12to2;
      gz+= K2*difu2.z + K1*difu.z + K0;
    } else {
      ftype efrom = e_fr[iold].z;
      ftype efrom2 = efrom*efrom;
      ftype eto = e_to[inew].z;
      ftype K0 =  efrom*(efrom2-e12to2+eto*eto)*(efrom+eto);
      ftype K1 =  (2*efrom+eto)*(e12to2-eto*eto)-efrom2*(4*efrom+3*eto);
      ftype K2 =  6*efrom2+3*efrom*eto-(e12to2-eto*eto);
      gz+= K2*difu2.z + K1*difu.z + K0;
    }*/
    //ftype difM4moment = u0.x+e_fr[iold]; difM4moment = difM4moment*difM4moment*difM4moment*difM4moment;
    //difM4moment*=0;
    //if(fabs(M4)>0.001) difM4moment*=-0.5;
    //else difM4moment*=0;
    //const ftype3 gxyz = K4*difu4 + K3*difu3 + K2*difu2 + K1*difu + K0;
    const ftype gxyz = gx*gy*gz;
    ftype fadd=0,hadd=0;
    if(Dim==3) fadd=wk*wl*wm*gxyz*f[iold];
    if(!isnan(fadd) && !isinf(fadd)) fnew+= fadd;
    
//     if(ix==123 && inew==0) printf("K0=%g %g %g | e_fr[iold]=%g %g %g (iold=%d inew=%d) T01 = %g %g   |  fnew=%g\n",K0.x,K0.y,K0.z,ef[iold].x,ef[iold].y,ef[iold].z,iold,inew, T0,T1,  fnew);

    #ifdef DDF
    if(Dim==1) hadd=wk*gx*h[iold];
    if(!isnan(hadd) && !isinf(hadd)) hnew+= hadd;
    #endif
  }
  #else
  if(e[inew].x==0) wk= -1./T1;
  if(e[inew].y==0) wl= -1./T1;
  if(e[inew].z==0) wm= -1./T1;
  for(int iold=0; iold<Qn; iold++) {
    ftype Ax_i = (u1.x-u0.x)*sqrt(Tbase)-e[iold].x*sqrt(fabs(T0))*sign(T0);
    ftype Ay_j = (u1.y-u0.y)*sqrt(Tbase)-e[iold].y*sqrt(fabs(T0))*sign(T0);
    ftype Az_p = (u1.z-u0.z)*sqrt(Tbase)-e[iold].z*sqrt(fabs(T0))*sign(T0);
    ftype Bx_ik = T1; if(e[inew].x!=0) Bx_ik= e[inew].x*sqrt(fabs(T1))*sign(T1)*Ax_i;
    ftype By_jl = T1; if(e[inew].y!=0) By_jl= e[inew].y*sqrt(fabs(T1))*sign(T1)*Ay_j;
    ftype Bz_pm = T1; if(e[inew].z!=0) Bz_pm= e[inew].z*sqrt(fabs(T1))*sign(T1)*Az_p;
    ftype gx=Ax_i*Ax_i-Bx_ik;
    ftype gy=Ay_j*Ay_j-By_jl;
    ftype gz=Az_p*Az_p-Bz_pm;
    ftype fadd=0,hadd=0;
    if(Dim==1) fadd=wk*gx*f[iold]; else 
    if(Dim==3) fadd=wk*wl*wm*gx*gy*gz*f[iold];
    if(!isnan(fadd) && !isinf(fadd)) fnew+= fadd;
    #ifdef DDF
    if(Dim==1) hadd=wk*gx*h[iold]; else 
    if(Dim==3) hadd=wk*wl*wm*gx*gy*gz*h[iold];
    if(!isnan(hadd) && !isinf(hadd)) hnew+= hadd;
    #endif
  }
  #endif
  return make_ftype2(fnew,hnew);
}
inline void __device__ Cell::mirrX(){
  uT.x=-uT.x;
  ftype fnew[Qn];
  for(int i=0; i<Qn; i++) {
    if(e_c[i].x==0) fnew[i]=f[i];
    else fnew[i]=f[reverseX[i]];
  }
  for(int i=0; i<Qn; i++) f[i]=fnew[i];
}
inline void __device__ Cell::mirrY(){
  uT.y=-uT.y;
  ftype fnew[Qn];
  for(int i=0; i<Qn; i++) {
    if(e_c[i].y==0) fnew[i]=f[i];
    else fnew[i]=f[reverseY[i]];
  }
  for(int i=0; i<Qn; i++) f[i]=fnew[i];
}
inline void __device__ Cell::mirrZ(){
  uT.z=-uT.z;
  ftype fnew[Qn];
  for(int i=0; i<Qn; i++) {
    if(e_c[i].z==0) fnew[i]=f[i];
    else fnew[i]=f[reverseZ[i]];
  }
  for(int i=0; i<Qn; i++) f[i]=fnew[i];
}

inline void __device__ Cell::calc_eq(ftype feq[Qn]) {
  ftype3 u3=make_ftype3(uT.x,uT.y,uT.z);
  const ftype mxw0 = ftype(1) - dot(u3,u3)*ftype(0.5)*dcs2;
  for(int i=0; i<Qn; i++) {
    ftype3 eidx = make_ftype3(e_c[i]);
    ftype eu =  dot(eidx,u3)*dcs2;
    ftype mxw   = mxw0 + eu + eu*eu*ftype(0.5);
    feq[i] = w_c[i]*rho*mxw;
  }
}

inline void __device__ Cell::collision(){
#if defined D1Q3 || defined D1Q5 || defined D1Q7 || defined D3Q125
#else
  rho=0; 
  for(int i=0; i<Qn; i++) rho+= f[i];
 //for(int i=0; i<Qn; i++) f[i] = w_c[i]*rho; return;
  ftype3 u3=make_ftype3(0,0,0); 
  for(int i=0; i<Qn; i++) u3+= f[i]*make_ftype3(e[i]); if(rho!=0) u3/=rho;
  uT.x=u3.x; uT.y=u3.y; uT.z=u3.z;
  register ftype feq[Qn];
  calc_eq(feq);
  register ftype _f[Qn];
  for(int i=0; i<Qn; i++) _f[i] = f[i]-pars.pc.dtau*(f[i]-feq[i]);
  for(int i=0; i<Qn; i++) f[i] = _f[reverseXYZ[i]];
}
__forceinline__ __device__ void Cell::fast_collision(){
  rho=0; 
  register ftype3 u3=make_ftype3(0,0,0); 
  for(int i=0; i<Qn; i++) rho+= f[i];
  u3.x = E_MATRIX_X(f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
  u3.y = E_MATRIX_Y(f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
  u3.z = E_MATRIX_Z(f[0],f[1],f[2],f[3],f[4],f[5],f[6],f[7],f[8],f[9],f[10],f[11],f[12],f[13],f[14],f[15],f[16],f[17],f[18]);
//   if(rho!=0) u3/=rho;
  if(rho!=0) { u3.x=__fdividef(u3.x,rho); u3.y=__fdividef(u3.y,rho); u3.z=__fdividef(u3.z,rho); }
  uT.x=u3.x; uT.y=u3.y; uT.z=u3.z;
  register ftype feq,frev,eu;
  const register ftype mxw0 = ftype(1) - dot(u3,u3)*ftype(0.5)*dcs2;
  const ftype dtau = pars.pc.dtau;
  feq = W0*rho*mxw0; f[0]+= (feq-f[0])*dtau;
  feq = W1*rho*( u3.x*(u3.x*dcs4/2 - dcs2) + mxw0 ); frev  = (feq-f[1] )*dtau+f[1] ;
  feq = W1*rho*( u3.x*(u3.x*dcs4/2 + dcs2) + mxw0 ); f[1]  = (feq-f[2] )*dtau+f[2] ; f[2]=frev;
  feq = W1*rho*( u3.y*(u3.y*dcs4/2 - dcs2) + mxw0 ); frev  = (feq-f[4] )*dtau+f[4] ;
  feq = W1*rho*( u3.y*(u3.y*dcs4/2 + dcs2) + mxw0 ); f[4]  = (feq-f[7] )*dtau+f[7] ; f[7]=frev;
  feq = W1*rho*( u3.z*(u3.z*dcs4/2 - dcs2) + mxw0 ); frev  = (feq-f[10])*dtau+f[10];
  feq = W1*rho*( u3.z*(u3.z*dcs4/2 + dcs2) + mxw0 ); f[10] = (feq-f[13])*dtau+f[13]; f[13]=frev;
  eu=u3.x+u3.y; feq = W2*rho*( eu*(eu*dcs4/2 - dcs2) + mxw0 ); frev  = (feq-f[3 ])*dtau+f[3] ;
  eu=u3.x+u3.y; feq = W2*rho*( eu*(eu*dcs4/2 + dcs2) + mxw0 ); f[3]  = (feq-f[8 ])*dtau+f[8] ; f[8]=frev;
  eu=u3.x-u3.y; feq = W2*rho*( eu*(eu*dcs4/2 + dcs2) + mxw0 ); frev  = (feq-f[5 ])*dtau+f[5] ;
  eu=u3.x-u3.y; feq = W2*rho*( eu*(eu*dcs4/2 - dcs2) + mxw0 ); f[5]  = (feq-f[6 ])*dtau+f[6] ; f[6]=frev;
  eu=u3.x+u3.z; feq = W2*rho*( eu*(eu*dcs4/2 - dcs2) + mxw0 ); frev  = (feq-f[9 ])*dtau+f[9] ;
  eu=u3.x+u3.z; feq = W2*rho*( eu*(eu*dcs4/2 + dcs2) + mxw0 ); f[9]  = (feq-f[14])*dtau+f[14]; f[14]=frev;
  eu=u3.x-u3.z; feq = W2*rho*( eu*(eu*dcs4/2 + dcs2) + mxw0 ); frev  = (feq-f[11])*dtau+f[11] ;
  eu=u3.x-u3.z; feq = W2*rho*( eu*(eu*dcs4/2 - dcs2) + mxw0 ); f[11] = (feq-f[12])*dtau+f[12]; f[12]=frev;
  eu=u3.y+u3.z; feq = W2*rho*( eu*(eu*dcs4/2 - dcs2) + mxw0 ); frev  = (feq-f[15])*dtau+f[15] ;
  eu=u3.y+u3.z; feq = W2*rho*( eu*(eu*dcs4/2 + dcs2) + mxw0 ); f[15] = (feq-f[18])*dtau+f[18]; f[18]=frev;
  eu=u3.y-u3.z; feq = W2*rho*( eu*(eu*dcs4/2 + dcs2) + mxw0 ); frev  = (feq-f[16])*dtau+f[16] ;
  eu=u3.y-u3.z; feq = W2*rho*( eu*(eu*dcs4/2 - dcs2) + mxw0 ); f[16] = (feq-f[17])*dtau+f[17]; f[17]=frev;
#endif
}

//----------------------------------------------------------------------------------//
inline static void __host__ __device__ setDFs(float f[Qn], float rho=1, ftype3 u=(ftype3){0,0,0}) {
  for(int i=0; i<Qn; i++) {
    ftype3 eidx = make_ftype3(e_c[i]);
    ftype eu =  dot(eidx,u)*dcs2;
    ftype u2=dot(u,u)*dcs2;
    ftype e2=dot(eidx,eidx)*dcs2;
    ftype mxw   = ( 1. + eu + eu*eu*0.5 - u2*0.5 );
    f[i] = w_c[i]*rho*mxw;
  }
}
inline static float4 __host__ __device__ getRV(float f[Qn]) {
  ftype r=0;ftype3 v=make_ftype3(0,0,0);
  for(int i=0; i<Qn; i++) {
    r+=f[i];
    v+=f[i]*make_ftype3(e_c[i]);
  }
  if(r!=0) v/=r;
  return make_float4(r,v.x,v.y,v.z);
}
inline static void __device__ fetchDFs(float f[Qn], int x, int y, int z) {
  for(int i=0; i<No; i++) {
    float2 fpair;
    //printf("%p\n",pars.data);
    surf3Dread(&fpair, pars.texSdata, (x*No+i)*sizeof(float2), y, z);
    f[2*i]=fpair.x;
    if(2*i+1<Qn) f[2*i+1]=fpair.y;
  }
}

int _main(int argc, char** argv);
