#include "structs.cuh"
#if 0
template<int Noreder> __device__ inline void LBMpullstep(const uint3 crd, const int it) {
  register Cell cell;
  cell.gather_group(x,y,z);
  cell.collision();
  const int2 gcind = Group::ind_conv(x, y, z);
  const int gindex=gcind.x, cindex=gcind.y;
  pars.groups[gindex]->pack(cell, cindex);
}

__global__ /*__launch_bounds__(brick.x*brick.y*brick.z)*/ void Layer_pull(){
  const uint3 bcrd = Sem_triple::unzip(blockIdx.x);
  const uint3 coord = threadIdx+brick*bcrd;
  LBMpullstep<2>(coord, pars.iStep);
}
#endif
