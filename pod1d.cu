#include "params.h"

const int Nx=512;
const int Nth=128;


__global__ __launch_bounds__(Nth) void LBMpull(){
  Cell c0; c0.load(coord.x,coord.y,coord.z);
  const ftype TbaseLat=TLat;
  if(c0.rho==0 && c0.uT.x==0 && c0.uT.y==0 && c0.uT.z==0 && c0.uT.w==0) c0.set_base(1,make_ftype4(0,0,0,TbaseLat));
  volatile int Niter=0;
  while(Niter<1000) {
    Niter++;
    ftype3 vi[Qn];
    const ftype4 uT=c0.uT;
    ftype3 ru = make_ftype3(uT.x,uT.y,uT.z);
    for(int i=0; i<Qn; i++) vi[i] = sqrt(uT.w/TLat)*make_ftype3(e[i])+ru;
    Cell cnew;
    for(int i=0; i<Qn; i++) {
      cnew.f[i] = propagate<2>(make_ftype3(coord)-vi[i], i, uT);
    }
    cnew.calcMoments(vi);
    if(cnew.uT.w<0) { 
//       if(coord.y==0 && coord.z==0) printf("neg T=%g prev=%g iter=%d\n",cnew.uT.w, c0.uT.w, Niter);
      cnew.uT.w=c0.uT.w;
    }
//     if(cnew.rho<0) cnew.rho=c0.rho;
//      if(coord.x==pars.Nx/2 && coord.y==pars.Ny/2 && coord.z==pars.Nz/2) 
//        printf("Niter=%d old_R_u_T=%g (%g %g %g) %g new=%g (%g %g %g) %g\n",Niter, c0.rho, c0.uT.x, c0.uT.y, c0.uT.z, c0.uT.w, cnew.rho, cnew.uT.x, cnew.uT.y, cnew.uT.z, cnew.uT.w);
//      if(coord.x==pars.Nx/2+1 && coord.y==pars.Ny/2 && coord.z==pars.Nz/2) 
//        printf("Niter=%d fi=%g %g %g new=%g %g %g\n",Niter, c0.f[8], c0.f[0], c0.f[1], cnew.f[8], cnew.f[0], cnew.f[1]);
//      if(coord.x==pars.Nx/2+1 && coord.y==pars.Ny/2 && coord.z==pars.Nz/2) 
//        printf("Niter=%d old_R_ux_T=%g (%g) %g new=%g (%g) %g\n",Niter, c0.rho, c0.uT.x, c0.uT.w, cnew.rho, cnew.uT.x, cnew.uT.w);
    pars.cinfo[Cinfo::ind_zip(coord.x,coord.y,coord.z)].niter=Niter;
    if(isConv(cnew,c0)) { c0=cnew; break; }
    c0=cnew;
  }
  atomicMax(pars.NiterMax,Niter);
  const ftype dtau = pars.pc.dtau;
  //const ftype dtau = 1./(0.2*TLat/c0.uT.w+0.5);
  ftype3 M1L = make_ftype3(0,0,0);
  ftype M2L = 0, M0L=0;
  ftype H=0;
  for(int i=0; i<Qn; i++) {
    M0L += c0.f[i];
    M1L += c0.f[i]*make_ftype3(e_c[i]*dx);
    M2L += c0.f[i]*dot(e_c[i]*dx,e_c[i]*dx);
    H+= c0.f[i]*log(c0.f[i]/w_c[i]);
  }
  if(coord.x==pars.Nx/2+1 && coord.y==pars.Ny/2 && coord.z==pars.Nz/2)
      printf("pre  collision R_ux_T=%g (%g) %g H=%g\n", M0L, M1L.x/M0L, (M2L-dot(M1L,M1L)/M0L)/Dim, H);
  for(int i=0; i<Qn; i++) {
    ftype feq = c0.rho*w[i];
    c0.f[i]+= dtau*(feq-c0.f[i]);
  }
  M1L = make_ftype3(0,0,0);
  M2L = 0; M0L=0; H=0;
  for(int i=0; i<Qn; i++) {
    M0L += c0.f[i];
    M1L += c0.f[i]*make_ftype3(e_c[i]*dx);
    M2L += c0.f[i]*dot(e_c[i]*dx,e_c[i]*dx);
    H+= c0.f[i]*log(c0.f[i]/w_c[i]);
  }
  if(coord.x==pars.Nx/2+1 && coord.y==pars.Ny/2 && coord.z==pars.Nz/2)
      printf("post-collision R_ux_T=%g (%g) %g H=%g\n", M0L, M1L.x/M0L, (M2L-dot(M1L,M1L)/M0L)/Dim, H);
  c0.save(coord.x,coord.y,coord.z);

}


void step(int it){

  LBMpull<<<Nx/Nth,Nth>>>(); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  
}

int main(int argc, char** argv) {
  assert(Nx%Nth);
  for(int it=0; it<Nt; it++) step(it);
  return 0;
}
