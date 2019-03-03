#include "structs.cuh"

#include "im2D.h"
#include "im3D.hpp"

#include "steps.cu"
#include "amr.cu"

template<int Norder> inline __device__ void LBMpull(const int3 coord);

__global__ __launch_bounds__(brick.x*brick.y*brick.z) void PonDstep(){
  const int3 coord = make_int3(blockIdx*brick+threadIdx);
  LBMpull<2>(coord);
}
__global__ __launch_bounds__(brick.x*brick.y*brick.z) void TimeCorr(ftype dtscale){
  const int3 coord = make_int3(blockIdx*brick+threadIdx);
  const int2 gcind = Group::ind_conv(coord.x, coord.y, coord.z);
  const int gindex=gcind.x, cindex=gcind.y;
  Cell c0;
  pars.storeGroups[gindex].unpack(c0, cindex);
  const ftype TbaseLat=TLat;
  c0.uT.x*= dtscale;
  c0.uT.y*= dtscale;
  c0.uT.z*= dtscale;
  c0.uT.w*= dtscale*dtscale;

  pars.storeGroups[gindex].pack(c0, cindex);
}

__global__ __launch_bounds__(brick.x*brick.y*brick.z) void cinfo_setnew_reset(){
   const int3 coord = make_int3(blockIdx*brick+threadIdx);
   Cinfo& cinf  = pars.cinfo[Cinfo::ind_zip(coord.x,coord.y,coord.z)];
   cinf.setnew = cinf.set;
}
__global__ __launch_bounds__(brick.x*brick.y*brick.z) void cinfo_set_new(){
   const int3 coord = make_int3(blockIdx*brick+threadIdx);
   Cinfo& cinf = pars.cinfo[Cinfo::ind_zip(coord.x,coord.y,coord.z)];
   cinf.set = cinf.setnew;
}
__global__ __launch_bounds__(brick.x*brick.y*brick.z) void UpdateMesh(){
   const int3 coord = make_int3(blockIdx*brick+threadIdx);
   Cinfo& cinf = pars.cinfo[Cinfo::ind_zip(coord.x,coord.y,coord.z)];
   Cell c0; c0.load(coord.x,coord.y,coord.z);
   Cell c1; c1.load_st(coord.x,coord.y,coord.z);
   if(cinf.set>0) {
     #ifdef NRMESH
     c0.p=c1.p;
     c0.save_ld(coord.x,coord.y,coord.z);
     #endif
   }
}
__global__ __launch_bounds__(brick.x*brick.y*brick.z) void calcTotMoments(){
   const int3 coord = make_int3(blockIdx*brick+threadIdx);
   Cinfo& cinf = pars.cinfo[Cinfo::ind_zip(coord.x,coord.y,coord.z)];
   Cell c0; c0.load_st(coord.x,coord.y,coord.z);
   Cell cm; cm.load_st((coord.x-1+pars.Nx)%pars.Nx,coord.y,coord.z);
   Cell cp; cp.load_st((coord.x+1+pars.Nx)%pars.Nx,coord.y,coord.z);
   if(cinf.set>0) {
     #ifdef NRMESH
     if(cp.p.x<cm.p.x) cp.p.x+=pars.Nx;
     const ftype dstep = 0.5*(cp.p.x-cm.p.x);
     #else 
     const ftype dstep = 1.0;
     #endif
     atomicAdd(pars.mass, c0.rho*dstep);
     atomicAdd((ftype*)pars.moment, c0.rho*c0.uT.x*dstep);
     atomicAdd(pars.enrg, c0.rho*(c0.uT.x*c0.uT.x+c0.uT.y*c0.uT.y+c0.uT.z*c0.uT.z+Dim*c0.uT.w)*dstep);
   }
}
__global__ __launch_bounds__(1) void print_data(int x){
   const int3 coord = make_int3(x,0,0);
   Cinfo& cinf = pars.cinfo[Cinfo::ind_zip(coord.x,coord.y,coord.z)];
   Cell c; c.load_st(coord.x,coord.y,coord.z);
   ftype cpnt=coord.x;
   #ifdef NRMESH
   cpnt = c.p.x;
   #endif
   if(cinf.set>0) printf("%d %.17f %.17f %.17f %.17f\n", coord.x, cpnt, c.rho, c.uT.x, c.uT.w);
}

extern im3D_pars im3DHost;

extern int maxR;
void drop(std::string* fprefix);
void diagnostics();
void calcStep(){
  //dim3 dimGrid((parsHost.Nx+dimBlock.x-1)/dimBlock.x, (parsHost.Ny+dimBlock.y-1)/dimBlock.y, (parsHost.Nz+dimBlock.z-1)/dimBlock.z);
    
  printf("Step %6d\n",parsHost.iStep);
  printf("running %d blocks with %d threads\n", tfNx*tfNy*tfNz, brick.x*brick.y*brick.z);
  cuTimer t0;
  CHECK_ERROR(cudaMemset(parsHost.NiterMax, 0, sizeof(int)));
  CHECK_ERROR(cudaMemset(parsHost.mass    , 0, sizeof(ftype)));
  CHECK_ERROR(cudaMemset(parsHost.enrg    , 0, sizeof(ftype)));
  CHECK_ERROR(cudaMemset(parsHost.moment  , 0, sizeof(ftype3)));
  uint3 block = dim3(parsHost.Nx/brick.x,parsHost.Ny/brick.y,parsHost.Nz/brick.z);
  #if defined D3Q125 && 1
  PonDstepD3Q125<<<dim3(parsHost.Nx,parsHost.Ny,parsHost.Nz),Qn>>>(); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  #else
  PonDstep  <<<block,brick>>>(); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  //UpdateMesh<<<block,brick>>>(); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  cinfo_setnew_reset<<<block,brick>>>(); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
   //AMR_add   <<<block,brick>>>(); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() ); 
  cinfo_set_new     <<<block,brick>>>(); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
   //AMR_remove<<<1,1>>>(); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() ); 
  #endif
  std::string drop_dir(im3DHost.drop_dir);
  drop(&drop_dir);
  diagnostics();

  parsHost.curtime+=parsHost.dt;
  ftype dtscale=1;
  
  if((parsHost.iStep+1)%1==0) {
    printf("\n\n#x rho ux T\n");
    for(int x=0; x<parsHost.Nx; x++) {
      //print_data<<<1,1>>>(x); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() ); 
    }
  }
//   usleep(10000);
//   if(*parsHost.NiterMax>=1000) dtscale=0.5;
//   else if(*parsHost.NiterMax<=4) dtscale=2;
//   if(dtscale!=1) {  TimeCorr<<<block,brick>>>(dtscale); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );  }
  CHECK_ERROR(cudaMemset(parsHost.mass    , 0, sizeof(ftype)));
  CHECK_ERROR(cudaMemset(parsHost.enrg    , 0, sizeof(ftype)));
  CHECK_ERROR(cudaMemset(parsHost.moment  , 0, sizeof(ftype3)));
  calcTotMoments<<<block,brick>>>(); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  parsHost.dt*=dtscale;
  double time1 = t0.gettime();
  printf("curtime=%g dt=%g  (dtscale %g, NiterMax=%d)\n", parsHost.curtime, parsHost.dt, dtscale, *parsHost.NiterMax);
  printf("time=%9g perf=%9g MLUps\n",time1, parsHost.Nx*parsHost.Ny*parsHost.Nz*Nt/time1*1e-3);
  printf("mass %.17f\n",*parsHost.mass);
  ftype3 moment = *parsHost.moment;
  printf("moment %.17f\n",moment.x);
  printf("enrg %.17f\n",*parsHost.enrg);


  parsHost.iStep++;
  parsHost.swapLS();
  //parsHost.reset();
  copy2dev( parsHost, pars );
}
