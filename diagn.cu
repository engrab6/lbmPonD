#include "structs.cuh"
__managed__ ftype k,eps;
__global__ void diagn(){
  // Calculate surface coordinates
  unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
  unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;
  if (x < pars.Nx && y < pars.Ny && z<pars.Nz) {
    Cell c; 
    const int2 gcind = Group::ind_conv(x, y, z);
    const int gindex=gcind.x, cindex=gcind.y;
    pars.storeGroups[gindex].unpack(c, cindex);
    ftype3 vel = make_ftype3(c.uT.x, c.uT.y, c.uT.z);
    ftype const d8pi3 = 1.0/(8*M_PI*M_PI*M_PI);
    const ftype dx = 2*M_PI/pars.Nx;
    const ftype dy = 2*M_PI/pars.Ny;
    const ftype dz = 2*M_PI/pars.Nz;
    atomicAdd(&k,0.5*d8pi3*(vel.x*vel.x+vel.y*vel.y+vel.z*vel.z)*dx*dy*dz);
    ftype visc=(1/pars.pc.dtau-0.5)*c.uT.w;
//    ftype epsadd = 2*visc*d8pi3*
  }
}

void diagnostics(){
  dim3 dimBlock(16, 16, 1);
  dim3 dimGrid((parsHost.Nx+dimBlock.x-1)/dimBlock.x, (parsHost.Ny+dimBlock.y-1)/dimBlock.y, (parsHost.Nz+dimBlock.z-1)/dimBlock.z);
  k=0;
  diagn<<<dimGrid, dimBlock>>>(); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() ); 
  printf("k= %g\n",k);
}
