#include "structs.cuh"
template<int odevity> __device__ __forceinline__ void  stepLBM(const int3 crd, const int it) {
//    return;
  int x=crd.x, y=crd.y, z=crd.z;
  if(x<0 || y<0 || z<0 || x>=pars.Nx || y>=pars.Ny || z>=pars.Nz) return;
  const long index = Cell::ind_zip(x,y,z);
  if(odevity==111110) {
    register Cell cell;
    register ftype cf[Qn];
//     cell.gather_group(x,y,z);
//     cell.fast_collision();
//     cell.scatter_group(x,y,z);
    int cx=x,cy=y,cz=z;
        #pragma unroll 19
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
          cf[i] = pars.loadGroups[gindex].f[inear][cindex];
        }
        #pragma unroll 19
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
          pars.storeGroups[gindex].f[inear][cindex] = cf[i];
        }
  } else {
    /*Cell* curcell = &pars.cells[index];
    Cell cellcopy = *curcell;
    //for(int i=0;i<Qn;i++) cellcopy.f[i]++;
    cellcopy.collision();
    *curcell = cellcopy;*/
    const int2 gcind = Group::ind_conv(x, y, z);
    const int gindex=gcind.x, cindex=gcind.y;
    Group* g = &pars.groups[gindex];
    register Cell cell; g->unpack(cell, cindex);
//     for(int i=0;i<Qn;i++) cell.f[i]++;
    cell.fast_collision();
    g->pack(cell, cindex);
  }
}

inline __device__ bool isConv(Cell& cnew, Cell& cold) {
  ftype val_new[] = {cnew.uT.x, cnew.uT.y, cnew.uT.z, cnew.uT.w};
  ftype val_old[] = {cold.uT.x, cold.uT.y, cold.uT.z, cold.uT.w};
  #ifdef USE_DOUBLE
  const ftype err_abs=1e-12;
  const ftype err_rel=1e-10;
  #else
  const ftype err_abs=1e-12;
  const ftype err_rel=1e-10;
  #endif
  for(int i=0; i<sizeof(val_new)/sizeof(val_new[0]); i++) {
    if(fabs(val_new[i]-val_old[i])>= err_abs+err_rel*fabs(val_new[i])) return false;
  }
  for(int i=0; i<Qn; i++) {
    if(fabs(cnew.f[i]-cold.f[i])>= err_abs+err_rel*fabs(cnew.f[i])) return false;
  }
  return true;
}
template<int N> inline __device__ ftype2 propagate     ( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir=0);
inline                 __device__ ftype2 propagateTVD2 ( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir=0);
template<int N> inline __device__ ftype2 propagateMacro( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir=0);
template<int N> inline __device__ ftype2 propagate_nonequal_1D( ftype3 crd, const int i, const ftype4 base_gauge, int niter=0 );
template<int N,int M> inline __device__ ftype2 propagateWENO( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir=0, const int niter=0);
template<int N> inline __device__ ftype2 propagateWENOmacro( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir=0, const int niter=0);
template<int N> inline __device__ ftype propagateTVD( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir=0);
inline __device__ ftype2 propagate3( ftype3 crd, int3 crd_left, const int i, const ftype4 base_gauge);
inline __device__ ftype2 propagateNR_1D(ftype3 crd, int3 crd_cnt,  const int i, const ftype4 base_gauge);
template<int N> inline __device__ ftype2 LagrEval( ftype3 crd, int3 crd_left, const int i, const ftype4 base_gauge);

inline ftype2 __device__ weno(ftype3 crdM2, ftype3 crdM1, ftype3 crd0, ftype3 crdP1, ftype3 crdP2, const int i, const ftype4 uT){
  ftype2 F1,F2,F3;
  ftype c1=1; 
//   if(crd0.x>crdP1.x) c1=-1;
  ftype2 Fm2 = c1*propagate<3>(crdM2, i, uT);
  ftype2 Fm1 = c1*propagate<3>(crdM1, i, uT);
  ftype2 F0  = c1*propagate<3>(crd0 , i, uT);
  ftype2 Fp1 = c1*propagate<3>(crdP1, i, uT);
  ftype2 Fp2 = c1*propagate<3>(crdP2, i, uT);

  F1 = 1./3.*Fm2 - 7./6.*Fm1 + 11./6.*F0;
  F2 = -1./6.*Fm1 + 5./6.*F0 + 1./3.*Fp1;
  F3 = 1./3.*F0 + 5./6.*Fp1 - 1./6.*Fp2;
  
  ftype delta1=1./10., delta2=3./5., delta3=3./10.;
  ftype2 sigma1 = 13./12.*(Fm2-2*Fm1+F0 )*(Fm2-2*Fm1+F0 ) + 1./4.*(Fm2-4*Fm1+3*F0)*(Fm2-4*Fm1+3*F0);
  ftype2 sigma2 = 13./12.*(Fm1-2*F0 +Fp1)*(Fm1-2*F0 +Fp1) + 1./4.*(Fm1-3*Fp1)*(Fm1-3*Fp1);
  ftype2 sigma3 = 13./12.*(F0 -2*Fp1+Fp2)*(F0 -2*Fp1+Fp2) + 1./4.*(3*F0-4*Fp1+Fp2)*(3*F0-4*Fp1+Fp2);
  ftype2 wq1; wq1.x = delta1/((1e-6+sigma1.x)*(1e-6+sigma1.x)); wq1.y = delta1/((1e-6+sigma1.y)*(1e-6+sigma1.y));
  ftype2 wq2; wq2.x = delta2/((1e-6+sigma2.x)*(1e-6+sigma2.x)); wq2.y = delta2/((1e-6+sigma2.y)*(1e-6+sigma2.y));
  ftype2 wq3; wq3.x = delta3/((1e-6+sigma3.x)*(1e-6+sigma3.x)); wq3.y = delta3/((1e-6+sigma3.y)*(1e-6+sigma3.y));
  ftype2 wsum = wq1+wq2+wq3;
  ftype2 w1 = wq1/wsum, w2 = wq2/wsum, w3 = wq3/wsum;
//   return Fp1-F0;
  return w1*F1+w2*F2+w3*F3;
}
inline  ftype __device__ minmod(const ftype a, const ftype b){
  if(a*b<=0) return 0; else
  if(fabs(a)<fabs(b)) return a;
  else return b;
}
inline ftype __device__ FluxTVDdif(ftype3 crdM2, ftype3 crdM1, ftype3 crd0, ftype3 crdP1, const int i, const ftype4 uT, ftype dir){
  ftype Fm1 = propagate<4>(crdM1 , i, uT).x;
  ftype Fm2 = propagate<4>(crdM2 , i, uT).x;
  ftype F0  = propagate<4>(crd0  , i, uT).x;
  ftype Fp1 = propagate<4>(crdP1 , i, uT).x;
  return F0-Fm1 + sign(dir)*0.5*minmod(F0-Fm1,Fp1-F0) - sign(dir)*0.5*minmod(F0-Fm1,Fm1-Fm2);
}
template<int N> inline __device__ ftype propagateTVD( ftype3 crd, const int i, const ftype4 uT, const ftype dir ){
  int3 pos = make_int3(floor(crd.x),floor(crd.x+0.5),floor(crd.x+0.5));
  if(dir>0) pos.x++;
  pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx; pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny; pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
  int3 pos0=pos, posM2=pos, posM1=pos, posP1=pos;
  posM2.x=pos0.x-2*sign(dir); posM2.x = (posM2.x%pars.Nx+pars.Nx)%pars.Nx;
  posM1.x=pos0.x-1*sign(dir); posM1.x = (posM1.x%pars.Nx+pars.Nx)%pars.Nx;
  posP1.x=pos0.x+1*sign(dir); posP1.x = (posP1.x%pars.Nx+pars.Nx)%pars.Nx;
  Cell c;
  c.load(pos0.x ,pos.y,pos.z); ftype F0  = c.transfer_gauge(uT, i).x;
  c.load(posM1.x,pos.y,pos.z); ftype Fm1 = c.transfer_gauge(uT, i).x;
  c.load(posM2.x,pos.y,pos.z); ftype Fm2 = c.transfer_gauge(uT, i).x;
  c.load(posP1.x,pos.y,pos.z); ftype Fp1 = c.transfer_gauge(uT, i).x;
  ftype difC = 0.5*minmod(Fp1-F0 ,F0-Fm1)*dir*0.01*0;
  ftype difU = 0.5*minmod(Fm1-Fm2,F0-Fm1)*dir*0.01*0;
//   if(pos.y==pars.Ny/2 && pos.z==pars.Nz/2 && i==0 && difC-difU!=0) printf("crd=%g, dir=%g, uT=(%g,%g) difC=%g difU=%g\n", crd.x, dir, uT.x, uT.w, difC,difU);
  return propagate<N>(crd, i, uT).x-difC+difU;
}
template<> __device__ ftype2 propagateWENO<4,4>( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir, const int);
template<> __device__ ftype2 propagateWENO<3,4>( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir, const int);
template<> __device__ ftype2 propagateWENO<3,5>( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir, const int);
template<> __device__ ftype2 propagateWENO<5,7>( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir, const int);
template<int N> inline __device__ ftype LagrPol(int ix,int iy,int iz, ftype sx, ftype sy,ftype sz);

const int N=3;
__shared__ Cell cfix[N][N][N];
__global__ __launch_bounds__(Qn) void PonDstepD3Q125(){
  const int3 coord=make_int3(blockIdx);
  const int inew=threadIdx.x;
  ftype dtscale=1;
  Cinfo& cinf = pars.cinfo[Cinfo::ind_zip(coord.x,coord.y,coord.z)];
  if(cinf.set==0) return;
  Cell c0; c0.load(coord.x,coord.y,coord.z);
  const ftype TbaseLat=TLat;
  if(c0.rho==0 && c0.uT.x==0 && c0.uT.y==0 && c0.uT.z==0 && c0.uT.w==0) c0.set_base(1,make_ftype4(0,0,0,TbaseLat));
  Cell cb=c0;
  volatile int Niter=0;
  const int N=3;
  const int3 crd_left = coord-make_int3(1,1,1);
  for(int icx=0; icx<N; icx++) for(int icy=0; icy<N; icy++) for(int icz=0; icz<N; icz++){
    int3 pos = crd_left+make_int3(icx,icy,icz);
    int doReverseX=0,doReverseY=0,doReverseZ=0;
    const bool outflow=0;
    if(outflow) {
      if(pos.x<0       ) { pos.x = -1-pos.x;                doReverseX=1; }
      if(pos.x>=pars.Nx) { pos.x = pars.Nx-1-pos.x+pars.Nx; doReverseX=1; }
      if(pos.y<0       ) { pos.y = -1-pos.y;                doReverseY=1; }
      if(pos.y>=pars.Ny) { pos.y = pars.Ny-1-pos.y+pars.Ny; doReverseY=1; }
      if(pos.z<0       ) { pos.z = -1-pos.z;                doReverseZ=1; }
      if(pos.z>=pars.Nz) { pos.z = pars.Nz-1-pos.z+pars.Nz; doReverseZ=1; }
    }
    pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx;
    pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny;
    pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
    const int2 gcind = Group::ind_conv(pos.x, pos.y, pos.z); const int gindex=gcind.x, cindex=gcind.y;
    Group* ldgrp = &pars.loadGroups[gindex];
    int iload=inew;
    if(doReverseX) iload=reverseX[iload];
    if(doReverseY) iload=reverseY[iload];
    if(doReverseZ) iload=reverseZ[iload];
    cfix[icx][icy][icz].f[inew]=ldgrp->f[iload][cindex];
    if(threadIdx.x==0) {
      cfix[icx][icy][icz].rho=ldgrp->rho[cindex];
      cfix[icx][icy][icz].uT =ldgrp->uT[cindex];
      if(doReverseX) cfix[icx][icy][icz].uT.x*=-1;
      if(doReverseY) cfix[icx][icy][icz].uT.y*=-1;
      if(doReverseZ) cfix[icx][icy][icz].uT.z*=-1;
    }
  }
//   if(coord.x==0 && coord.y==0 && coord.z==0) for(int i=0; i<Qn; i++) printf("reverse: %d : %d\n",i, reverseZ[i]);
  __syncthreads();
  while(Niter<100) {
    Niter++;
    ftype3 vi[Qn];
    const ftype4 uT=c0.uT;
    ftype3 ru = make_ftype3(uT.x,uT.y,uT.z);
    for(int i=0; i<Qn; i++) vi[i] = sqrt(uT.w/TLat)*ef[i]+ru;
    __shared__ Cell cnew;
    ftype2 newval[N][N][N];
    for(int icx=0; icx<N; icx++) for(int icy=0; icy<N; icy++) for(int icz=0; icz<N; icz++) {
      newval[icx][icy][icz] = cfix[icx][icy][icz].transfer_gauge(uT, inew, (icx==1 && icy==1 && icz==1 && coord.x==pars.Nx/2 && coord.y==pars.Ny/2 && coord.z==pars.Nz/2)?123:0);
    }
    const ftype3 crd = make_ftype3(coord)-vi[inew];
    const ftype3 shift = crd-make_ftype3(crd_left);
    ftype2 fpr=make_ftype2(0,0);
    for(int ix=0; ix<N; ix++) for(int iy=0; iy<N; iy++) for(int iz=0; iz<N; iz++) {
      fpr+= LagrPol<N>(ix,iy,iz,shift.x,shift.y,shift.z)*newval[ix][iy][iz];
    }
    ftype2 newfh = fpr;
    cnew.f[inew] = newfh.x;
    #ifdef DDF
    cnew.h[inew] = newfh.y;
    #endif
    __syncthreads();
    if(threadIdx.x==0) {
      cnew.calcMoments(vi,uT.x);
      if(cnew.uT.w<0) cnew.uT.w=c0.uT.w;
      cinf.niter=Niter;
    }
    __syncthreads();
    if(isConv(cnew,c0)) { c0=cnew; break; }
    c0=cnew;
    __syncthreads();
  }
  atomicMax(pars.NiterMax,Niter);
  ftype dtau = pars.pc.dtau;

  ftype feq = c0.rho*w[inew];
  ftype heq = c0.rho*w[inew]*(Dim*c0.uT.w+c0.uT.x*c0.uT.x+c0.uT.y*c0.uT.y+c0.uT.z*c0.uT.z);
  c0.f[inew]+= dtau*(feq-c0.f[inew]);
  #ifdef DDF
  c0.h[imew]+= dtau*(heq-c0.h[inew]);
  #endif
  
  const int2 gcind = Group::ind_conv(coord.x, coord.y, coord.z); const int gindex=gcind.x, cindex=gcind.y;
  Group* stgrp = &pars.storeGroups[gindex];
  stgrp->f[inew][cindex]=c0.f[inew];
  if(threadIdx.x==0) {
    stgrp->rho[cindex]=c0.rho;
    stgrp->uT[cindex]=c0.uT;
    atomicAdd(pars.mass, c0.rho);
    atomicAdd((ftype*)pars.moment, c0.rho*c0.uT.x);
    atomicAdd(pars.enrg, c0.rho*(c0.uT.x*c0.uT.x+c0.uT.y*c0.uT.y+c0.uT.z*c0.uT.z+Dim*c0.uT.w));
  }
  __syncthreads();
}
__global__ __launch_bounds__(Qn) void PonDstepD3Q125_iterative(){
  const int3 coord=make_int3(blockIdx);
  const int inew=threadIdx.x;
  ftype dtscale=1;
  Cinfo& cinf = pars.cinfo[Cinfo::ind_zip(coord.x,coord.y,coord.z)];
  if(cinf.set==0) return;
  Cell c0; c0.load(coord.x,coord.y,coord.z);
  const ftype TbaseLat=TLat;
  if(c0.rho==0 && c0.uT.x==0 && c0.uT.y==0 && c0.uT.z==0 && c0.uT.w==0) c0.set_base(1,make_ftype4(0,0,0,TbaseLat));
  Cell cb=c0;
  volatile int Niter=0;
  const int N=3;
  const int3 crd_left = coord-make_int3(1,1,1);
  for(int icx=0; icx<N; icx++) for(int icy=0; icy<N; icy++) for(int icz=0; icz<N; icz++){
    int3 pos = crd_left+make_int3(icx,icy,icz);
    int doReverseX=0,doReverseY=0,doReverseZ=0;
    const bool outflow=0;
    if(outflow) {
      if(pos.x<0       ) { pos.x = -1-pos.x;                doReverseX=1; }
      if(pos.x>=pars.Nx) { pos.x = pars.Nx-1-pos.x+pars.Nx; doReverseX=1; }
      if(pos.y<0       ) { pos.y = -1-pos.y;                doReverseY=1; }
      if(pos.y>=pars.Ny) { pos.y = pars.Ny-1-pos.y+pars.Ny; doReverseY=1; }
      if(pos.z<0       ) { pos.z = -1-pos.z;                doReverseZ=1; }
      if(pos.z>=pars.Nz) { pos.z = pars.Nz-1-pos.z+pars.Nz; doReverseZ=1; }
    }
    pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx;
    pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny;
    pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
    const int2 gcind = Group::ind_conv(pos.x, pos.y, pos.z); const int gindex=gcind.x, cindex=gcind.y;
    Group* ldgrp = &pars.loadGroups[gindex];
    cfix[icx][icy][icz].f[inew]=ldgrp->f[inew][cindex];
    if(threadIdx.x==0) {
      cfix[icx][icy][icz].rho=ldgrp->rho[cindex];
      cfix[icx][icy][icz].uT =ldgrp->uT[cindex];
      if(doReverseX) cfix[icx][icy][icz].mirrX();
      if(doReverseY) cfix[icx][icy][icz].mirrY();
      if(doReverseZ) cfix[icx][icy][icz].mirrZ();
    }
  }
  __syncthreads();
  while(Niter<100) {
    Niter++;
    ftype3 vi[Qn];
    const ftype4 uT=c0.uT;
    ftype3 ru = make_ftype3(uT.x,uT.y,uT.z);
    for(int i=0; i<Qn; i++) vi[i] = sqrt(uT.w/TLat)*ef[i]+ru;
    __shared__ Cell cnew;
    ftype2 newval[N][N][N];
    for(int icx=0; icx<N; icx++) for(int icy=0; icy<N; icy++) for(int icz=0; icz<N; icz++) {
      newval[icx][icy][icz] = cfix[icx][icy][icz].transfer_gauge(uT, inew, (icx==1 && icy==1 && icz==1 && coord.x==pars.Nx/2 && coord.y==pars.Ny/2 && coord.z==pars.Nz/2)?123:0);
    }
    const ftype3 crd = make_ftype3(coord)-vi[inew];
    const ftype3 shift = crd-make_ftype3(crd_left);
    ftype2 fpr=make_ftype2(0,0);
    for(int ix=0; ix<N; ix++) for(int iy=0; iy<N; iy++) for(int iz=0; iz<N; iz++) {
      fpr+= LagrPol<N>(ix,iy,iz,shift.x,shift.y,shift.z)*newval[ix][iy][iz];
    }
    ftype2 newfh = fpr;
    cnew.f[inew] = newfh.x;
    #ifdef DDF
    cnew.h[inew] = newfh.y;
    #endif
    __syncthreads();
    if(threadIdx.x==0) {
      cnew.calcMoments(vi,uT.x);
      if(cnew.uT.w<0) cnew.uT.w=c0.uT.w;
      cinf.niter=Niter;
    }
    __syncthreads();
    if(isConv(cnew,c0)) { c0=cnew; break; }
    c0=cnew;
    __syncthreads();
  }
  atomicMax(pars.NiterMax,Niter);
  ftype dtau = pars.pc.dtau;

  ftype feq = c0.rho*w[inew];
  ftype heq = c0.rho*w[inew]*(Dim*c0.uT.w+c0.uT.x*c0.uT.x+c0.uT.y*c0.uT.y+c0.uT.z*c0.uT.z);
  c0.f[inew]+= dtau*(feq-c0.f[inew]);
  #ifdef DDF
  c0.h[imew]+= dtau*(heq-c0.h[inew]);
  #endif
  
  const int2 gcind = Group::ind_conv(coord.x, coord.y, coord.z); const int gindex=gcind.x, cindex=gcind.y;
  Group* stgrp = &pars.storeGroups[gindex];
  stgrp->f[inew][cindex]=c0.f[inew];
  if(threadIdx.x==0) {
    stgrp->rho[cindex]=c0.rho;
    stgrp->uT[cindex]=c0.uT;
    atomicAdd(pars.mass, c0.rho);
    atomicAdd((ftype*)pars.moment, c0.rho*c0.uT.x);
    atomicAdd(pars.enrg, c0.rho*(c0.uT.x*c0.uT.x+c0.uT.y*c0.uT.y+c0.uT.z*c0.uT.z+Dim*c0.uT.w));
  }
  __syncthreads();
//   if(coord.x==0 && coord.y==0 && coord.z==0) {
//     Group* ldgrp = &pars.loadGroups[gindex];
//     ftype divu=0;
//     for(int inda=0; inda<3; inda++) for(int indb=0; indb<3; indb++) divu+= 0.5*(cfix[2][ia][ib].uT.x-cfix[0][ia][ib].uT.x);
//     for(int inda=0; inda<3; inda++) for(int indb=0; indb<3; indb++) divu+= 0.5*(cfix[ia][2][ib].uT.y-cfix[ia][0][ib].uT.y);
//     for(int inda=0; inda<3; inda++) for(int indb=0; indb<3; indb++) divu+= 0.5*(cfix[ia][ib][2].uT.z-cfix[ia][ib][0].uT.z);
//     for(int inda=0; inda<3; inda++) for(int indb=0; indb<3; indb++) divu+= 0.5*(cfix[2][ia][ib].uT.x-cfix[0][ia][ib].uT.x);
//     ftype assum_rho= ldgrp->rho[cindex]-0.5*(cfix[2][0][0].uT.x-cfix[0][0][0].uT.x);
//     printf("rnew=%g assumed=%g\n", stgrp->rho[cindex], assum_rho);
//   }
}
template<int Norder> __forceinline__ __device__ void LBMpull(const int3 coord){
  ftype dtscale=1;
  Cinfo& cinf = pars.cinfo[Cinfo::ind_zip(coord.x,coord.y,coord.z)];
  if(cinf.set==0) return;
  Cell c0; c0.load(coord.x,coord.y,coord.z);
  Cell c0_st; c0_st.load_st(coord.x,coord.y,coord.z);
  const ftype TbaseLat=TLat;
  if(c0.rho==0 && c0.uT.x==0 && c0.uT.y==0 && c0.uT.z==0 && c0.uT.w==0) c0.set_base(1,make_ftype4(0,0,0,TbaseLat));
  Cell cb=c0;
  volatile int Niter=0;
  #if defined D3Q125 && 1
  const int N=3;
  const int3 crd_left = coord-make_int3(1,1,1);
  Cell cfix[N][N][N];
  for(int icx=0; icx<N; icx++) for(int icy=0; icy<N; icy++) for(int icz=0; icz<N; icz++){
    int3 pos = crd_left+make_int3(icx,icy,icz);
    int doReverse=0;
    const bool outflow=0;
    if(outflow) {
      if(pos.x<0       ) pos.x = -1-pos.x;
      if(pos.x>=pars.Nx) pos.x = pars.Nx-1-pos.x+pars.Nx;
    }
    pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx;
    pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny;
    pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
    cfix[icx][icy][icz].load(pos.x,pos.y,pos.z);
//     if(doReverse) cfix[icx][icy][icz].reverse();
  }
  while(Niter<100) {
    Niter++;
    ftype3 vi[Qn];
    const ftype4 uT=c0.uT;
    ftype3 ru = make_ftype3(uT.x,uT.y,uT.z);
    for(int i=0; i<Qn; i++) vi[i] = sqrt(uT.w/TLat)*ef[i]+ru;
//     Cell c[N][N][N];
//     for(int icx=0; icx<N; icx++) for(int icy=0; icy<N; icy++) for(int icz=0; icz<N; icz++) c[icx][icy][icz]=cfix[icx][icy][icz];
    Cell cnew;
    #ifdef NRMESH
    cnew.p=c0.p;
    #endif
    ftype2 newval[Qn][N][N][N];
    for(int icx=0; icx<N; icx++) for(int icy=0; icy<N; icy++) for(int icz=0; icz<N; icz++) {
      for(int i=0; i<Qn; i++) newval[i][icx][icy][icz] = cfix[icx][icy][icz].transfer_gauge(uT, i, (icx==1 && icy==1 && icz==1 && coord.x==pars.Nx/2 && coord.y==pars.Ny/2 && coord.z==pars.Nz/2)?123:0);
    }
    for(int i=0; i<Qn; i++) {
      const ftype3 crd = make_ftype3(coord)-vi[i];
      const ftype3 shift = crd-make_ftype3(crd_left);
      ftype2 fpr=make_ftype2(0,0);
      for(int ix=0; ix<N; ix++) for(int iy=0; iy<N; iy++) for(int iz=0; iz<N; iz++) {
        fpr+= LagrPol<N>(ix,iy,iz,shift.x,shift.y,shift.z)*newval[i][ix][iy][iz];
      }
      ftype2 newfh = fpr;
      cnew.f[i] = newfh.x;
      #ifdef DDF
      cnew.h[i] = newfh.y;
      #endif
//       if(coord.x==pars.Nx/2 && coord.y==pars.Ny/2 && coord.z==pars.Nz/2) printf("iter=%d i=%d newval=%g\n", Niter, i, newval[i][1][1][1].x);
    }
    cnew.calcMoments(vi,uT.x);
//     if(coord.x==pars.Nx/2 && coord.y==pars.Ny/2 && coord.z==pars.Nz/2) printf("iter=%d c0: %g (%g %g %g) %g;  cnew=%g (%g %g %g) %g\n", Niter, c0.rho,c0.uT.x,c0.uT.y,c0.uT.z,c0.uT.w, cnew.rho,cnew.uT.x,cnew.uT.y,cnew.uT.z,cnew.uT.w);
    if(cnew.uT.w<0) cnew.uT.w=c0.uT.w;
    cinf.niter=Niter;
    if(isConv(cnew,c0)) { c0=cnew; break; }
    c0=cnew;
  }
  #else
  while(Niter<100) {
    Niter++;
    ftype3 vi[Qn];
    const ftype4 uT=c0.uT;
    ftype3 ru = make_ftype3(uT.x,uT.y,uT.z);
    #if defined D1Q5 || defined D3Q125
    for(int i=0; i<Qn; i++) vi[i] = sqrt(uT.w/TLat)*ef[i]+ru;
    #else
    for(int i=0; i<Qn; i++) vi[i] = sqrt(uT.w/TLat)*make_ftype3(e[i])+ru;
    #endif
    Cell cnew;
    #ifdef NRMESH
    cnew.p=c0.p;
    #endif
    ftype M4=0;
    Cell cm,cp; 
    cm.load((coord.x-1+pars.Nx)%pars.Nx,coord.y,coord.z);
    cp.load((coord.x+1+pars.Nx)%pars.Nx,coord.y,coord.z);
    M4 = (cp.rho+cm.rho-2*c0.rho)/c0.rho;
    //for(int i=0; i<Qn; i++) M4+=c0.f[i]*vi[i].x*vi[i].x*vi[i].x*vi[i].x;
    for(int i=0; i<Qn; i++) {
      ftype2 newfh=make_ftype2(0,0);
//       cnew.f[i] = propagateWENO<5>(make_ftype3(coord)-vi[i], i, uT, vi[i].x);
//       newfh = propagate<3>(make_ftype3(coord), i, uT) - ( weno(make_ftype3(coord)-vi[i]*2, make_ftype3(coord)-vi[i]*1, make_ftype3(coord)        , make_ftype3(coord)+vi[i]*1, make_ftype3(coord)+vi[i]*2, i, uT) -
//                                                           weno(make_ftype3(coord)-vi[i]*3, make_ftype3(coord)-vi[i]*2, make_ftype3(coord)-vi[i]*1, make_ftype3(coord)        , make_ftype3(coord)+vi[i]*1, i, uT)  );
//       cnew.f[i] = propagate<4>(make_ftype3(coord), i, uT) - ( FluxTVDdif(make_ftype3(coord)-vi[i]*2, make_ftype3(coord)-vi[i]*1, make_ftype3(coord)+vi[i]*0, make_ftype3(coord)+vi[i]*1, i, uT, vi[i].x));
//      cnew.f[i] = propagateTVD<10>(make_ftype3(coord)-vi[i], i, uT, vi[i].x);
//       newfh = propagate_nonequal_1D<2>(make_ftype3(coord)-vi[i], i, uT);
//        newfh = propagate     <3>(make_ftype3(coord)-vi[i], i, uT,-vi[i].x);
       newfh = propagate     <4>(make_ftype3(coord)-vi[i], i, uT,M4);
//       newfh = propagate3(make_ftype3(coord)-vi[i],coord-make_int3(1,1,1), i, uT);
//        newfh = propagateNR_1D(c0_st.p-vi[i],coord, i, uT);
//       newfh = propagateTVD2(make_ftype3(coord)-vi[i], i, uT,M4);
//       newfh = propagateMacro<4>(make_ftype3(coord)-vi[i], i, uT,-vi[i].x);

//        newfh = propagateWENO<5>(make_ftype3(coord)-vi[i], i, uT, vi[i].x, Niter);
//          newfh = propagateWENO<3,5>(make_ftype3(coord)-vi[i], i, uT, vi[i].x, Niter);
//           newfh = propagateWENO<3,4>(make_ftype3(coord)-vi[i], i, uT, vi[i].x, Niter);
//       newfh = propagateWENO<5,7>(make_ftype3(coord)-vi[i], i, uT, vi[i].x, Niter);
//          newfh = propagateWENO<4,4>(make_ftype3(coord)-vi[i], i, uT, vi[i].x, Niter);

//       newfh = propagateWENOmacro<4>(make_ftype3(coord)-vi[i], i, uT, vi[i].x, Niter);
      cnew.f[i] = newfh.x;
      #ifdef DDF
      cnew.h[i] = newfh.y;
      #endif
    }
    cnew.calcMoments(vi,uT.x);
/*ftype4 newGauge=cnew.uT; ftype newRho=cnew.rho;
cnew.uT=c0.uT;
for(int i=0; i<Qn; i++) cnew.f[i] = cnew.transfer_gauge(cnew.uT, i);
c0.uT=newGauge; c0.rho=newRho;
break;*/
    if(cnew.uT.w<0) { 
//       if(coord.y==0 && coord.z==0) printf("neg T=%g prev=%g iter=%d\n",cnew.uT.w, c0.uT.w, Niter);
      cnew.uT.w=c0.uT.w;
    }
//     if(cnew.rho<0) cnew.rho=c0.rho;
      if(coord.x==pars.Nx/2 && coord.y==pars.Ny/2 && coord.z==pars.Nz/2) { 
        Cell ct; ct.load(coord.x+1,coord.y,coord.z);
        /*printf("Niter=%d old_R_u_T=%g (%g) %g new=%g (%g) %g | f0=/%g, %g, %g/ | fnew=/%g, %g, %g/ |\n        cr= %g (%g) %g | frght=/%g,%g,%g/->/%g,%g,%g/\n",
        Niter, c0.rho, c0.uT.x, c0.uT.w, cnew.rho, cnew.uT.x, cnew.uT.w, c0.f[2],c0.f[0],c0.f[1], cnew.f[2],cnew.f[0],cnew.f[1],
        ct.rho, ct.uT.x, ct.uT.w,
        ct.f[2],ct.f[0],ct.f[1], ct.transfer_gauge(c0.uT, 2),ct.transfer_gauge(c0.uT, 0),ct.transfer_gauge(c0.uT, 1)
        );*/
        /*ftype ctM0=ct.f[0]+ct.f[1]+ct.f[2]+ct.f[3]+ct.f[4];
        ftype ctM1=ct.f[0]*(ct.uT.x)+ct.f[1]*(ec1*sqrt(ct.uT.w/TLat)+ct.uT.x)+ct.f[2]*(-ec1*sqrt(ct.uT.w/TLat)+ct.uT.x)+ct.f[3]*(ec2*sqrt(ct.uT.w/TLat)+ct.uT.x)+ct.f[4]*(-ec2*sqrt(ct.uT.w/TLat)+ct.uT.x);
        ftype ctM2=ct.f[0]*(ct.uT.x)*(ct.uT.x)+ct.f[1]*(ec1*sqrt(ct.uT.w/TLat)+ct.uT.x)*(ec1*sqrt(ct.uT.w/TLat)+ct.uT.x)+ct.f[2]*(-ec1*sqrt(ct.uT.w/TLat)+ct.uT.x)*(-ec1*sqrt(ct.uT.w/TLat)+ct.uT.x)+ct.f[3]*(ec2*sqrt(ct.uT.w/TLat)+ct.uT.x)*(ec2*sqrt(ct.uT.w/TLat)+ct.uT.x)+ct.f[4]*(-ec2*sqrt(ct.uT.w/TLat)+ct.uT.x)*(-ec2*sqrt(ct.uT.w/TLat)+ct.uT.x);
        ftype ctM3=ct.f[0]*(ct.uT.x)*(ct.uT.x)*(ct.uT.x)+ct.f[1]*(ec1*sqrt(ct.uT.w/TLat)+ct.uT.x)*(ec1*sqrt(ct.uT.w/TLat)+ct.uT.x)*(ec1*sqrt(ct.uT.w/TLat)+ct.uT.x)+ct.f[2]*(-ec1*sqrt(ct.uT.w/TLat)+ct.uT.x)*(-ec1*sqrt(ct.uT.w/TLat)+ct.uT.x)*(-ec1*sqrt(ct.uT.w/TLat)+ct.uT.x)+ct.f[3]*(ec2*sqrt(ct.uT.w/TLat)+ct.uT.x)*(ec2*sqrt(ct.uT.w/TLat)+ct.uT.x)*(ec2*sqrt(ct.uT.w/TLat)+ct.uT.x)+ct.f[4]*(-ec2*sqrt(ct.uT.w/TLat)+ct.uT.x)*(-ec2*sqrt(ct.uT.w/TLat)+ct.uT.x)*(-ec2*sqrt(ct.uT.w/TLat)+ct.uT.x);
        ftype newf0 = ct.transfer_gauge(c0.uT, 0);
        ftype newf1 = ct.transfer_gauge(c0.uT, 1);
        ftype newf2 = ct.transfer_gauge(c0.uT, 2);
        ftype newf3 = ct.transfer_gauge(c0.uT, 3);
        ftype newf4 = ct.transfer_gauge(c0.uT, 4);
        ftype newM0 = newf0+newf1+newf2+newf3+newf4; 
        ftype newM1 = newf0*(c0.uT.x) + newf1*(ec1*sqrt(c0.uT.w/TLat)+c0.uT.x) + newf2*(-ec1*sqrt(c0.uT.w/TLat)+c0.uT.x) + newf3*(ec2*sqrt(c0.uT.w/TLat)+c0.uT.x) + newf4*(-ec2*sqrt(c0.uT.w/TLat)+c0.uT.x); 
        ftype newM2 = newf0*(c0.uT.x)*(c0.uT.x) + newf1*(ec1*sqrt(c0.uT.w/TLat)+c0.uT.x)*(ec1*sqrt(c0.uT.w/TLat)+c0.uT.x) + newf2*(-ec1*sqrt(c0.uT.w/TLat)+c0.uT.x)*(-ec1*sqrt(c0.uT.w/TLat)+c0.uT.x) + newf3*(ec2*sqrt(c0.uT.w/TLat)+c0.uT.x)*(ec2*sqrt(c0.uT.w/TLat)+c0.uT.x) + newf4*(-ec2*sqrt(c0.uT.w/TLat)+c0.uT.x)*(-ec2*sqrt(c0.uT.w/TLat)+c0.uT.x); 
        ftype newM3 = newf0*(c0.uT.x)*(c0.uT.x)*(c0.uT.x) + newf1*(ec1*sqrt(c0.uT.w/TLat)+c0.uT.x)*(ec1*sqrt(c0.uT.w/TLat)+c0.uT.x)*(ec1*sqrt(c0.uT.w/TLat)+c0.uT.x) + newf2*(-ec1*sqrt(c0.uT.w/TLat)+c0.uT.x)*(-ec1*sqrt(c0.uT.w/TLat)+c0.uT.x)*(-ec1*sqrt(c0.uT.w/TLat)+c0.uT.x) + newf3*(ec2*sqrt(c0.uT.w/TLat)+c0.uT.x)*(ec2*sqrt(c0.uT.w/TLat)+c0.uT.x)*(ec2*sqrt(c0.uT.w/TLat)+c0.uT.x) + newf4*(-ec2*sqrt(c0.uT.w/TLat)+c0.uT.x)*(-ec2*sqrt(c0.uT.w/TLat)+c0.uT.x)*(-ec2*sqrt(c0.uT.w/TLat)+c0.uT.x); 
        printf("Niter=%d ctM0=%g ctM1=%g ctM2=%g ctM3=%g newM0=%g newM1=%g newM2=%g newM3=%g f0=%.18f f1=%.18f f2=%.18f f3=%.18f f4=%.18f from=%.18f %.18f to=%.18f %.18f\n ",
        Niter,ctM0,ctM1,ctM2,ctM3,newM0,newM1,newM2,newM3,
        ct.f[0],ct.f[1],ct.f[2],ct.f[3],ct.f[4],ct.uT.x,ct.uT.w,c0.uT.x,c0.uT.w);*/
      }
//      if(abs(coord.x-1025)<1.5 && coord.y==pars.Ny/2 && coord.z==pars.Nz/2) 
//        printf("Niter=%d fi=%g %g %g new=%g %g %g\n",Niter, c0.f[8], c0.f[0], c0.f[1], cnew.f[8], cnew.f[0], cnew.f[1]);
//      if(abs(coord.x-1024)<1.5 && coord.y==pars.Ny/2 && coord.z==pars.Nz/2) 
//        printf("Niter=%d old_R_ux_T=%g (%g) %g new=%g (%g) %g\n",Niter, c0.rho, c0.uT.x, c0.uT.w, cnew.rho, cnew.uT.x, cnew.uT.w);
    cinf.niter=Niter;
    if(isConv(cnew,c0)) { c0=cnew; break; }
    c0=cnew;
  }
  #endif
  atomicMax(pars.NiterMax,Niter);
  ftype dtau = pars.pc.dtau;
  //const ftype dtau = 1./(0.2*TLat/c0.uT.w+0.5);
  Cell cm,cp;
  cm.load((coord.x-1+pars.Nx)%pars.Nx,coord.y,coord.z);
  cp.load((coord.x+1+pars.Nx)%pars.Nx,coord.y,coord.z);
  //if(fabs((cp.rho-cm.rho)/c0.rho)>0.1) dtau=1./200.;
  //if((cp.rho-c0.rho)*(c0.rho-cm.rho)<0) dtau=1./2.;
  ftype3 M1L = make_ftype3(0,0,0);
  ftype M2L = 0, M0L=0;
  ftype H=0;
  for(int i=0; i<Qn; i++) {
    M0L += c0.f[i];
    M1L += c0.f[i]*make_ftype3(e_c[i]*dx);
    M2L += c0.f[i]*dot(e_c[i]*dx,e_c[i]*dx);
    H+= c0.f[i]*log(c0.f[i]/w[i]);
  }
  if(coord.x==pars.Nx/2+1 && coord.y==pars.Ny/2 && coord.z==pars.Nz/2)
      printf("pre  collision R_ux_T=%g (%g) %g H=%g\n", M0L, M1L.x/M0L, (M2L-dot(M1L,M1L)/M0L)/Dim, H);
  const ftype cT=TLat;
  //ftype M4c = ec1d*ec1d*(c0.f[1]+c0.f[2])+ec2d*ec2d*(c0.f[3]+c0.f[4]);
  //M4c/=c0.rho;
  //M4c = 3*TLat*TLat-1*(M4c-3*TLat*TLat);
  //if((c0.f[1]-c0.f[2])*(c0.f[3]-c0.f[4])<0) printf("Vyvert f1-1 and f2-2: cell %d f=%.17f %.17f %.17f %.17f %.17f u,T=%.17f,%.17f\n",coord.x, c0.f[0],c0.f[1],c0.f[2],c0.f[3],c0.f[4], c0.uT.x, c0.uT.w);
  /*const ftype W0tmp=(M4c-(ec1d+ec2d)*cT+ec1d*ec2d)/(ec1d*ec2d);
  const ftype W1tmp=(cT*ec2d-M4c)/(2*ec1d*(ec2d-ec1d));
  const ftype W2tmp=(cT*ec1d-M4c)/(2*ec2d*(ec1d-ec2d));
  const ftype wTmp[Qn] = {W0tmp,W1tmp,W1tmp,W2tmp,W2tmp};*/
  ftype M3 = 0,M4=0;

  /*ftype vi[Qn];  for(int i=0; i<Qn; i++) vi[i] = sqrt(c0.uT.w/TLat)*ef[i].x+c0.uT.x;
  for(int i=0;i<Qn;i++) M3+= c0.f[i]*vi[i]*vi[i]*vi[i];
  for(int i=0;i<Qn;i++) M4+= c0.f[i]*vi[i]*vi[i]*vi[i]*vi[i];
  M4 = c0.rho*(3*TLat*TLat+6*c0.uT.x*c0.uT.x+c0.uT.x*c0.uT.x*c0.uT.x*c0.uT.x);
  const ftype A = M3/c0.rho-c0.uT.x*(3*c0.uT.w+c0.uT.x*c0.uT.x);
  const ftype B = (M4-4*c0.uT.x*M3)/c0.rho+3*c0.uT.x*c0.uT.x*(2*c0.uT.w+c0.uT.x*c0.uT.x);
  const ftype Tdrel = TLat/cT;
  const ftype W0tmp=(B*Tdrel*Tdrel-(ec1d+ec2d)*cT+ec1d*ec2d)/(ec1d*ec2d);
  const ftype W1Ptmp=(cT*ec2d-B*Tdrel*Tdrel-A*ec1*sqrt(Tdrel*Tdrel*Tdrel))/(2*ec1d*(ec2d-ec1d));
  const ftype W1Mtmp=(cT*ec2d-B*Tdrel*Tdrel+A*ec1*sqrt(Tdrel*Tdrel*Tdrel))/(2*ec1d*(ec2d-ec1d));
  const ftype W2Ptmp=(cT*ec1d-B*Tdrel*Tdrel-A*ec2*sqrt(Tdrel*Tdrel*Tdrel))/(2*ec2d*(ec1d-ec2d));
  const ftype W2Mtmp=(cT*ec1d-B*Tdrel*Tdrel+A*ec2*sqrt(Tdrel*Tdrel*Tdrel))/(2*ec2d*(ec1d-ec2d));
  const ftype wTmp[Qn] = {W0tmp,W1Ptmp,W1Mtmp,W2Ptmp,W2Mtmp};*/

  for(int i=0; i<Qn; i++) {
    ftype feq = c0.rho*w[i];
    ftype heq = c0.rho*w[i]*(Dim*c0.uT.w+c0.uT.x*c0.uT.x+c0.uT.y*c0.uT.y+c0.uT.z*c0.uT.z);
    c0.f[i]+= dtau*(feq-c0.f[i]);
    #ifdef DDF
    c0.h[i]+= dtau*(heq-c0.h[i]);
    #endif
  }
  #ifdef NRMESH
  //c0.p=c0_st.p+make_ftype3(c0.uT.x,c0.uT.y,c0.uT.z);
  //if(c0.p.x<0); c0.p.x+=pars.Nx;
  //if(c0.p.x>=pars.Nx); c0.p.x-=pars.Nx;
  #endif

  M1L = make_ftype3(0,0,0);
  ftype Hold=H;
  M2L = 0; M0L=0; H=0;
  for(int i=0; i<Qn; i++) {
    M0L += c0.f[i];
    M1L += c0.f[i]*make_ftype3(e_c[i]*dx);
    M2L += c0.f[i]*dot(e_c[i]*dx,e_c[i]*dx);
    H+= c0.f[i]*log(c0.f[i]/w[i]);
  }
  if(coord.x==pars.Nx/2+1 && coord.y==pars.Ny/2 && coord.z==pars.Nz/2)
      printf("post-collision R_ux_T=%g (%g) %g H=%g\n", M0L, M1L.x/M0L, (M2L-dot(M1L,M1L)/M0L)/Dim, H);
  c0.save(coord.x,coord.y,coord.z);
  atomicAdd(pars.mass, c0.rho);
  atomicAdd((ftype*)pars.moment, c0.rho*c0.uT.x);
  atomicAdd(pars.enrg, c0.rho*(c0.uT.x*c0.uT.x+c0.uT.y*c0.uT.y+c0.uT.z*c0.uT.z+Dim*c0.uT.w));
//  pars.cinfo[Cinfo::ind_zip(coord.x,coord.y,coord.z)].niter=H;
//   if(Hold<H) printf("ix=%d  H-functions: %.17f %.17f\n",coord.x, Hold, H);
}

template<int N> inline __device__ ftype LagrPolNR_1D(int ix, ftype sx, ftype bp[N]){
  ftype a=1;
  for(int ixp=0; ixp<N; ixp++) if(ixp!=ix) a*= (sx-bp[ixp])/(bp[ix]-bp[ixp]);
  return a;
}
template<int N> inline __device__ ftype LagrPol(int ix,int iy,int iz, ftype sx, ftype sy,ftype sz){
  ftype a=1;
  for(int ixp=0; ixp<N; ixp++) if(ixp!=ix) a*= (sx-ixp)/(ix-ixp);
  for(int iyp=0; iyp<N; iyp++) if(iyp!=iy) a*= (sy-iyp)/(iy-iyp);
  for(int izp=0; izp<N; izp++) if(izp!=iz) a*= (sz-izp)/(iz-izp);
  return a;
}
inline ftype __device__ LinPol(int ix,int iy,int iz, ftype sx, ftype sy,ftype sz){
  ftype a=1;
  if(ix==0) a*= 1-sx; else if(ix==1) a*=sx;
  if(iy==0) a*= 1-sy; else if(iy==1) a*=sy;
  if(iz==0) a*= 1-sz; else if(iz==1) a*=sz;
  return a;
}

template<int N, int M> inline __device__ ftype2 propagateWENO( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir, const int niter){
  int3 crd_left;
  if(N!=5) return make_ftype2(0,0);
  ftype lshift=0;
  //if(dir<0) lshift=1;
  lshift=1;
  ftype2 S0 = propagate3(crd, make_int3(floor(crd.x+0.5)-lshift+1, floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype2 S1 = propagate3(crd, make_int3(floor(crd.x+0.5)-lshift  , floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype2 S2 = propagate3(crd, make_int3(floor(crd.x+0.5)-lshift-1, floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype coeff3L=1, coeff3R=1, coeff5L=1, coeff5R=1;
  ftype cdist = crd.x-(int(floor(crd.x+0.5)));
  ftype d0=(cdist*cdist+3*cdist+2)/12.;
  ftype d2=(cdist*cdist-3*cdist+2)/12.;
  ftype d1=1.-d0-d2;
  crd_left = make_int3(floor(crd.x+0.5)-lshift-1,floor(crd.x+0.5),floor(crd.x+0.5));
  int3 pos = crd_left+make_int3(0,0,0);
  int3 pos1 = crd_left+make_int3(0,0,0);
  int3 pos2 = crd_left+make_int3(1,0,0);
  int3 pos3 = crd_left+make_int3(2,0,0);
  int3 pos4 = crd_left+make_int3(3,0,0);
  int3 pos5 = crd_left+make_int3(4,0,0);
  pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx; pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny; pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
  pos1.x= (pos1.x%pars.Nx+pars.Nx)%pars.Nx; pos1.y= (pos1.y%pars.Ny+pars.Ny)%pars.Ny; pos1.z = (pos1.z%pars.Nz+pars.Nz)%pars.Nz;
  pos2.x= (pos2.x%pars.Nx+pars.Nx)%pars.Nx; pos2.y= (pos2.y%pars.Ny+pars.Ny)%pars.Ny; pos2.z = (pos2.z%pars.Nz+pars.Nz)%pars.Nz;
  pos3.x= (pos3.x%pars.Nx+pars.Nx)%pars.Nx; pos3.y= (pos3.y%pars.Ny+pars.Ny)%pars.Ny; pos3.z = (pos3.z%pars.Nz+pars.Nz)%pars.Nz;
  pos4.x= (pos4.x%pars.Nx+pars.Nx)%pars.Nx; pos4.y= (pos4.y%pars.Ny+pars.Ny)%pars.Ny; pos4.z = (pos4.z%pars.Nz+pars.Nz)%pars.Nz;
  pos5.x= (pos5.x%pars.Nx+pars.Nx)%pars.Nx; pos5.y= (pos5.y%pars.Ny+pars.Ny)%pars.Ny; pos5.z = (pos5.z%pars.Nz+pars.Nz)%pars.Nz;
  Cell c1,c2,c3,c4,c5;
  c1.load(pos1.x,pos.y,pos.z); ftype u1 = c1.transfer_gauge(base_gauge, i).x;
  c2.load(pos2.x,pos.y,pos.z); ftype u2 = c2.transfer_gauge(base_gauge, i).x;
  c3.load(pos3.x,pos.y,pos.z); ftype u3 = c3.transfer_gauge(base_gauge, i).x;
  c4.load(pos4.x,pos.y,pos.z); ftype u4 = c4.transfer_gauge(base_gauge, i).x;
  c5.load(pos5.x,pos.y,pos.z); ftype u5 = c5.transfer_gauge(base_gauge, i).x;

  ftype beta0 = (10*u3*u3-31*u3*u4+25*u4*u4+11*u3*u5-19*u4*u5+4*u5*u5)/3.0;
  ftype beta1 = (4*u2*u2 -13*u2*u3+13*u3*u3+5*u2*u4 -13*u3*u4+4*u4*u4)/3.0;
  ftype beta2 = (4*u1*u1 -19*u1*u2+25*u2*u2+11*u1*u3-31*u2*u3+10*u3*u3)/3.0;
//   if(pars.iStep<100) { beta0=0; beta1=0; beta2=0; }
  ftype w0 = d0/((1e-6+beta0)*(1e-6+beta0));
  ftype w1 = d1/((1e-6+beta1)*(1e-6+beta1));
  ftype w2 = d2/((1e-6+beta2)*(1e-6+beta2));
  ftype wsum = w0+w1+w2;
  w0/= wsum; w1/= wsum; w2/= wsum;
  return w0*S0+w1*S1+w2*S2;
}
template<> inline __device__ ftype2 propagateWENO<3,4>( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir, const int niter){
  int3 crd_left;
  ftype2 S0 = LagrEval<3>(crd, make_int3(floor(crd.x)-1, floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype2 S1 = LagrEval<3>(crd, make_int3(floor(crd.x)  , floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype2 S01= LagrEval<2>(crd, make_int3(floor(crd.x), floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype cdist = crd.x-(int(floor(crd.x)));
  ftype d0=(2-cdist)/3.;
  ftype d1=1-d0; //(cdist+1)/3.;
  
  crd_left = make_int3(floor(crd.x)-1,floor(crd.x+0.5),floor(crd.x+0.5));
  int3 pos = crd_left+make_int3(0,0,0);
  int3 pos1 = crd_left+make_int3(0,0,0);
  int3 pos2 = crd_left+make_int3(1,0,0);
  int3 pos3 = crd_left+make_int3(2,0,0);
  int3 pos4 = crd_left+make_int3(3,0,0);
  pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx; pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny; pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
  pos1.x= (pos1.x%pars.Nx+pars.Nx)%pars.Nx; pos1.y= (pos1.y%pars.Ny+pars.Ny)%pars.Ny; pos1.z = (pos1.z%pars.Nz+pars.Nz)%pars.Nz;
  pos2.x= (pos2.x%pars.Nx+pars.Nx)%pars.Nx; pos2.y= (pos2.y%pars.Ny+pars.Ny)%pars.Ny; pos2.z = (pos2.z%pars.Nz+pars.Nz)%pars.Nz;
  pos3.x= (pos3.x%pars.Nx+pars.Nx)%pars.Nx; pos3.y= (pos3.y%pars.Ny+pars.Ny)%pars.Ny; pos3.z = (pos3.z%pars.Nz+pars.Nz)%pars.Nz;
  pos4.x= (pos4.x%pars.Nx+pars.Nx)%pars.Nx; pos4.y= (pos4.y%pars.Ny+pars.Ny)%pars.Ny; pos4.z = (pos4.z%pars.Nz+pars.Nz)%pars.Nz;
  Cell c1,c2,c3,c4;
  c1.load(pos1.x,pos.y,pos.z); ftype2 u1 = c1.transfer_gauge(base_gauge, i);
  c2.load(pos2.x,pos.y,pos.z); ftype2 u2 = c2.transfer_gauge(base_gauge, i);
  c3.load(pos3.x,pos.y,pos.z); ftype2 u3 = c3.transfer_gauge(base_gauge, i);
  c4.load(pos4.x,pos.y,pos.z); ftype2 u4 = c4.transfer_gauge(base_gauge, i);

  ftype2 beta0 = (25*u3*u3+(26*u1-76*u2)*u3+64*u2*u2-52*u1*u2+13*u1*u1)/12.0;
  ftype2 beta1 = (13*u4*u4+(26*u2-52*u3)*u4+64*u3*u3-76*u2*u3+25*u2*u2)/12.0;
  
//   if(floor(crd.x)==floor(crd.x+0.5)) {
//     beta0 = (4*u3*u3+(5*u1-13*u2)*u3+13*u2*u2-13*u1*u2+4*u1*u1)/3.0;
//     beta1 = (4*u4*u4+(11*u2-19*u3)*u4+25*u3*u3-31*u2*u3+10*u2*u2)/3.0;
//   } else {
//     beta0 = (10*u3*u3+(11*u1-31*u2)*u3+25*u2*u2-19*u1*u2+4*u1*u1)/3.0;
//     beta1 = (4*u4*u4+(5*u2-13*u3)*u4+13*u3*u3-13*u2*u3+4*u2*u2)/3.0;
//   }

//   if(pars.iStep<100) { beta0=0; beta1=0; beta2=0; }
  ftype2 w0 = d0*make_ftype2(1,1);
  ftype2 w1 = d1*make_ftype2(1,1);
  
  ftype smooth_coeff0=1e-2;
  ftype smooth_coeff1=1e-2;
  //if(floor(crd.x)==floor(crd.x+0.5)) smooth_coeff1=1e+1; else smooth_coeff0=1e+1;

  w0.x/= (smooth_coeff0+beta0.x)*(smooth_coeff0+beta0.x); w0.y/= (smooth_coeff0+beta0.y)*(smooth_coeff0+beta0.y);
  w1.x/= (smooth_coeff1+beta1.x)*(smooth_coeff1+beta1.x); w1.y/= (smooth_coeff1+beta1.y)*(smooth_coeff1+beta1.y);
  
  //if(fabs(beta0.x)>fabs(beta1.x)) {w0.x=0; w1.x=1;} else {w1.x=0; w0.x=1;}
  /*if(floor(crd.x)==floor(crd.x+0.5)) {
    if(fabs(beta0.x)>2*fabs(beta1.x)) {w0.x=0; w1.x=1; S1=S01; } else {w0.x=1; w1.x=0;}
  } else {
    if(fabs(beta1.x)>2*fabs(beta0.x)) {w0.x=0; w1.x=1; S1=S01; } else {w0.x=0; w1.x=1;}
  }*/
  //if(floor(crd.x)==floor(crd.x+0.5)   && i==0) { w0.x=0; w1.x=1; }
  //if(floor(crd.x)==floor(crd.x+0.5)-1 && i==0) { w0.x=1; w1.x=0; }
     
  ftype2 wsum = w0+w1;
  w0/= wsum; w1/= wsum;

  //if(int(crd.y)==0 && int(crd.z==0) && fabs(crd.x+dir-1024)<=4 && i==0)
  //printf("i=%d niter=%d crd.x=%.10f, w012=%g %g %g, d012=%g,%g,%g, beta012=%g,%g,%g u123456=%g,%g,%g,%g,%g,%g S012=%g %g %g\n",i,niter,crd.x, w0,w1,w2, d0,d1,d2, beta0,beta1,beta2, u1,u2,u3,u4,u5,u6, S0,S1,S2);
  //printf("i=%d niter=%d crd.x=%.10f, w012=%g %g %g, S012=%g %g %g\n",i,niter,crd.x, w0,w1,w2, S0,S1,S2);

  return w0*S0+w1*S1;
}
template<> inline __device__ ftype2 propagateWENO<5,7>( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir, const int niter){
  int3 crd_left;
  ftype2 S0 = LagrEval<5>(crd, make_int3(floor(crd.x+0.5)-3, floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype2 S1 = LagrEval<5>(crd, make_int3(floor(crd.x+0.5)-2, floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype2 S2 = LagrEval<5>(crd, make_int3(floor(crd.x+0.5)-1, floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype cdist = crd.x-(int(floor(crd.x+0.5)));
  ftype d0=(cdist*cdist-5*cdist+6)/30.;
  ftype d2=(cdist*cdist+5*cdist+6)/30.;
  ftype d1=1.-d0-d2;
  crd_left = make_int3(floor(crd.x+0.5)-3,floor(crd.x+0.5),floor(crd.x+0.5));
  int3 pos = crd_left+make_int3(0,0,0);
  int3 pos1 = crd_left+make_int3(0,0,0);
  int3 pos2 = crd_left+make_int3(1,0,0);
  int3 pos3 = crd_left+make_int3(2,0,0);
  int3 pos4 = crd_left+make_int3(3,0,0);
  int3 pos5 = crd_left+make_int3(4,0,0);
  int3 pos6 = crd_left+make_int3(5,0,0);
  int3 pos7 = crd_left+make_int3(6,0,0);
  pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx; pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny; pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
  pos1.x= (pos1.x%pars.Nx+pars.Nx)%pars.Nx; pos1.y= (pos1.y%pars.Ny+pars.Ny)%pars.Ny; pos1.z = (pos1.z%pars.Nz+pars.Nz)%pars.Nz;
  pos2.x= (pos2.x%pars.Nx+pars.Nx)%pars.Nx; pos2.y= (pos2.y%pars.Ny+pars.Ny)%pars.Ny; pos2.z = (pos2.z%pars.Nz+pars.Nz)%pars.Nz;
  pos3.x= (pos3.x%pars.Nx+pars.Nx)%pars.Nx; pos3.y= (pos3.y%pars.Ny+pars.Ny)%pars.Ny; pos3.z = (pos3.z%pars.Nz+pars.Nz)%pars.Nz;
  pos4.x= (pos4.x%pars.Nx+pars.Nx)%pars.Nx; pos4.y= (pos4.y%pars.Ny+pars.Ny)%pars.Ny; pos4.z = (pos4.z%pars.Nz+pars.Nz)%pars.Nz;
  pos5.x= (pos5.x%pars.Nx+pars.Nx)%pars.Nx; pos5.y= (pos5.y%pars.Ny+pars.Ny)%pars.Ny; pos5.z = (pos5.z%pars.Nz+pars.Nz)%pars.Nz;
  pos6.x= (pos6.x%pars.Nx+pars.Nx)%pars.Nx; pos6.y= (pos6.y%pars.Ny+pars.Ny)%pars.Ny; pos6.z = (pos6.z%pars.Nz+pars.Nz)%pars.Nz;
  pos7.x= (pos7.x%pars.Nx+pars.Nx)%pars.Nx; pos7.y= (pos7.y%pars.Ny+pars.Ny)%pars.Ny; pos7.z = (pos7.z%pars.Nz+pars.Nz)%pars.Nz;
  Cell c1,c2,c3,c4,c5,c6,c7;
  c1.load(pos1.x,pos.y,pos.z); ftype2 u1 = c1.transfer_gauge(base_gauge, i);
  c2.load(pos2.x,pos.y,pos.z); ftype2 u2 = c2.transfer_gauge(base_gauge, i);
  c3.load(pos3.x,pos.y,pos.z); ftype2 u3 = c3.transfer_gauge(base_gauge, i);
  c4.load(pos4.x,pos.y,pos.z); ftype2 u4 = c4.transfer_gauge(base_gauge, i);
  c5.load(pos5.x,pos.y,pos.z); ftype2 u5 = c5.transfer_gauge(base_gauge, i);
  c6.load(pos6.x,pos.y,pos.z); ftype2 u6 = c6.transfer_gauge(base_gauge, i);
  c7.load(pos7.x,pos.y,pos.z); ftype2 u7 = c7.transfer_gauge(base_gauge, i);

  ftype2 beta0 = (218654*u5*u5+(-1230721*u4+1288227*u3-595723*u2+100909*u1)*u5+1964729*u4*u4+(-4454616*u3+2119342*u2-363463*u1)*u4+2676879*u3*u3+(460407*u1-2647776*u2)*u3+682889*u2*u2-241621*u1*u2+21884*u1*u1)/60480.0;
  ftype2 beta1 = (21884*u6*u6+(-117931*u5+74217*u4+22727*u3-22781*u2)*u6+261209*u5*u5+(-592716*u4+165502*u3+22727*u2)*u5+518499*u4*u4+(74217*u2-592716*u3)*u4+261209*u3*u3-117931*u2*u3+21884*u2*u2)/60480.0;
  ftype2 beta2 = (21884*u7*u7+(-241621*u6+460407*u5-363463*u4+100909*u3)*u7+682889*u6*u6+(-2647776*u5+2119342*u4-595723*u3)*u6+2676879*u5*u5+(1288227*u3-4454616*u4)*u5+1964729*u4*u4-1230721*u3*u4+218654*u3*u3)/60480.0;

  ftype2 w0 = d0*make_ftype2(1,1);
  ftype2 w1 = d1*make_ftype2(1,1);
  ftype2 w2 = d2*make_ftype2(1,1);
  
  ftype smooth_coeff=1e-6;
  w0.x/= (smooth_coeff+beta0.x)*(smooth_coeff+beta0.x); w0.y/= (smooth_coeff+beta0.y)*(smooth_coeff+beta0.y);
  w1.x/= (smooth_coeff+beta1.x)*(smooth_coeff+beta1.x); w1.y/= (smooth_coeff+beta1.y)*(smooth_coeff+beta1.y);
  w2.x/= (smooth_coeff+beta2.x)*(smooth_coeff+beta2.x); w2.y/= (smooth_coeff+beta2.y)*(smooth_coeff+beta2.y);
  
  ftype2 wsum = w0+w1+w2;
  w0/= wsum; w1/= wsum; w2/= wsum;

  return w0*S0+w1*S1+w2*S2;
}
template<> inline __device__ ftype2 propagateWENO<3,5>( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir, const int niter){
  int3 crd_left;
  ftype2 S0 = LagrEval<3>(crd, make_int3(floor(crd.x+0.5)-2, floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype2 S1 = LagrEval<3>(crd, make_int3(floor(crd.x+0.5)-1, floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype2 S2 = LagrEval<3>(crd, make_int3(floor(crd.x+0.5)  , floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype cdist = crd.x-(int(floor(crd.x+0.5)));
  ftype d0=(cdist*cdist-3*cdist+2)/12.;
  ftype d2=(cdist*cdist+3*cdist+2)/12.;
  ftype d1=1.-d0-d2;
  crd_left = make_int3(floor(crd.x+0.5)-2,floor(crd.x+0.5),floor(crd.x+0.5));
  int3 pos = crd_left+make_int3(0,0,0);
  int3 pos1 = crd_left+make_int3(0,0,0);
  int3 pos2 = crd_left+make_int3(1,0,0);
  int3 pos3 = crd_left+make_int3(2,0,0);
  int3 pos4 = crd_left+make_int3(3,0,0);
  int3 pos5 = crd_left+make_int3(4,0,0);
  pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx; pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny; pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
  pos1.x= (pos1.x%pars.Nx+pars.Nx)%pars.Nx; pos1.y= (pos1.y%pars.Ny+pars.Ny)%pars.Ny; pos1.z = (pos1.z%pars.Nz+pars.Nz)%pars.Nz;
  pos2.x= (pos2.x%pars.Nx+pars.Nx)%pars.Nx; pos2.y= (pos2.y%pars.Ny+pars.Ny)%pars.Ny; pos2.z = (pos2.z%pars.Nz+pars.Nz)%pars.Nz;
  pos3.x= (pos3.x%pars.Nx+pars.Nx)%pars.Nx; pos3.y= (pos3.y%pars.Ny+pars.Ny)%pars.Ny; pos3.z = (pos3.z%pars.Nz+pars.Nz)%pars.Nz;
  pos4.x= (pos4.x%pars.Nx+pars.Nx)%pars.Nx; pos4.y= (pos4.y%pars.Ny+pars.Ny)%pars.Ny; pos4.z = (pos4.z%pars.Nz+pars.Nz)%pars.Nz;
  pos5.x= (pos5.x%pars.Nx+pars.Nx)%pars.Nx; pos5.y= (pos5.y%pars.Ny+pars.Ny)%pars.Ny; pos5.z = (pos5.z%pars.Nz+pars.Nz)%pars.Nz;
  Cell c1,c2,c3,c4,c5;
  c1.load(pos1.x,pos.y,pos.z); ftype2 u1 = c1.transfer_gauge(base_gauge, i);
  c2.load(pos2.x,pos.y,pos.z); ftype2 u2 = c2.transfer_gauge(base_gauge, i);
  c3.load(pos3.x,pos.y,pos.z); ftype2 u3 = c3.transfer_gauge(base_gauge, i);
  c4.load(pos4.x,pos.y,pos.z); ftype2 u4 = c4.transfer_gauge(base_gauge, i);
  c5.load(pos5.x,pos.y,pos.z); ftype2 u5 = c5.transfer_gauge(base_gauge, i);

  ftype2 beta0 = (10*u3*u3+(11*u1-31*u2)*u3+25*u2*u2-19*u1*u2+4*u1*u1)/3.0;
  ftype2 beta1 = (4*u4*u4+(5*u2-13*u3)*u4+13*u3*u3-13*u2*u3+4*u2*u2)/3.0;
  ftype2 beta2 = (4*u5*u5+(11*u3-19*u4)*u5+25*u4*u4-31*u3*u4+10*u3*u3)/3.0;

//   if(pars.iStep<100) { beta0=0; beta1=0; beta2=0; }
  ftype2 w0 = d0*make_ftype2(1,1);
  ftype2 w1 = d1*make_ftype2(1,1);
  ftype2 w2 = d2*make_ftype2(1,1);
  
  ftype smooth_coeff=1e-6;
  w0.x/= (smooth_coeff+beta0.x)*(smooth_coeff+beta0.x); w0.y/= (smooth_coeff+beta0.y)*(smooth_coeff+beta0.y);
  w1.x/= (smooth_coeff+beta1.x)*(smooth_coeff+beta1.x); w1.y/= (smooth_coeff+beta1.y)*(smooth_coeff+beta1.y);
  w2.x/= (smooth_coeff+beta2.x)*(smooth_coeff+beta2.x); w2.y/= (smooth_coeff+beta2.y)*(smooth_coeff+beta2.y);
  
  ftype2 wsum = w0+w1+w2;
  w0/= wsum; w1/= wsum; w2/= wsum;

  //if(int(crd.y)==0 && int(crd.z==0) && fabs(crd.x+dir-1024)<=4 && i==0)
  //printf("i=%d niter=%d crd.x=%.10f, w012=%g %g %g, d012=%g,%g,%g, beta012=%g,%g,%g u123456=%g,%g,%g,%g,%g,%g S012=%g %g %g\n",i,niter,crd.x, w0,w1,w2, d0,d1,d2, beta0,beta1,beta2, u1,u2,u3,u4,u5,u6, S0,S1,S2);
  //printf("i=%d niter=%d crd.x=%.10f, w012=%g %g %g, S012=%g %g %g\n",i,niter,crd.x, w0,w1,w2, S0,S1,S2);

  return w0*S0+w1*S1+w2*S2;
}

template<> inline __device__ ftype2 propagateWENO<4,4>( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir, const int niter){
  int3 crd_left;
  ftype2 S0 = LagrEval<4>(crd, make_int3(floor(crd.x)  , floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype2 S1 = LagrEval<4>(crd, make_int3(floor(crd.x)-1, floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype2 S2 = LagrEval<4>(crd, make_int3(floor(crd.x)-2, floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype cdist = crd.x-(int(floor(crd.x)));
  ftype d0=(cdist*cdist+3*cdist+2)/20.;
  ftype d2=(cdist*cdist-5*cdist+6)/20.;
  ftype d1=1.-d0-d2;
  crd_left = make_int3(floor(crd.x)-2,floor(crd.x+0.5),floor(crd.x+0.5));
  int3 pos = crd_left+make_int3(0,0,0);
  int3 pos1 = crd_left+make_int3(0,0,0);
  int3 pos2 = crd_left+make_int3(1,0,0);
  int3 pos3 = crd_left+make_int3(2,0,0);
  int3 pos4 = crd_left+make_int3(3,0,0);
  int3 pos5 = crd_left+make_int3(4,0,0);
  int3 pos6 = crd_left+make_int3(5,0,0);
  pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx; pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny; pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
  pos1.x= (pos1.x%pars.Nx+pars.Nx)%pars.Nx; pos1.y= (pos1.y%pars.Ny+pars.Ny)%pars.Ny; pos1.z = (pos1.z%pars.Nz+pars.Nz)%pars.Nz;
  pos2.x= (pos2.x%pars.Nx+pars.Nx)%pars.Nx; pos2.y= (pos2.y%pars.Ny+pars.Ny)%pars.Ny; pos2.z = (pos2.z%pars.Nz+pars.Nz)%pars.Nz;
  pos3.x= (pos3.x%pars.Nx+pars.Nx)%pars.Nx; pos3.y= (pos3.y%pars.Ny+pars.Ny)%pars.Ny; pos3.z = (pos3.z%pars.Nz+pars.Nz)%pars.Nz;
  pos4.x= (pos4.x%pars.Nx+pars.Nx)%pars.Nx; pos4.y= (pos4.y%pars.Ny+pars.Ny)%pars.Ny; pos4.z = (pos4.z%pars.Nz+pars.Nz)%pars.Nz;
  pos5.x= (pos5.x%pars.Nx+pars.Nx)%pars.Nx; pos5.y= (pos5.y%pars.Ny+pars.Ny)%pars.Ny; pos5.z = (pos5.z%pars.Nz+pars.Nz)%pars.Nz;
  pos6.x= (pos6.x%pars.Nx+pars.Nx)%pars.Nx; pos6.y= (pos6.y%pars.Ny+pars.Ny)%pars.Ny; pos6.z = (pos6.z%pars.Nz+pars.Nz)%pars.Nz;
  Cell c1,c2,c3,c4,c5,c6;
  c1.load(pos1.x,pos.y,pos.z); ftype2 u1 = c1.transfer_gauge(base_gauge, i);
  c2.load(pos2.x,pos.y,pos.z); ftype2 u2 = c2.transfer_gauge(base_gauge, i);
  c3.load(pos3.x,pos.y,pos.z); ftype2 u3 = c3.transfer_gauge(base_gauge, i);
  c4.load(pos4.x,pos.y,pos.z); ftype2 u4 = c4.transfer_gauge(base_gauge, i);
  c5.load(pos5.x,pos.y,pos.z); ftype2 u5 = c5.transfer_gauge(base_gauge, i);
  c6.load(pos6.x,pos.y,pos.z); ftype2 u6 = c6.transfer_gauge(base_gauge, i);

  ftype2 beta0 = (64*u6*u6+(-579*u5+774*u4-323*u3)*u6+1356*u5*u5+(1554*u3-3687*u4)*u5+2706*u4*u4-2499*u3*u4+634*u3*u3)/180.0;
  ftype2 beta1 = (64*u5*u5+(-189*u4-6*u3+67*u2)*u5+366*u4*u4+(-537*u3-6*u2)*u4+366*u3*u3-189*u2*u3+64*u2*u2)/180.0;
  ftype2 beta2 = (634*u4*u4+(-2499*u3+1554*u2-323*u1)*u4+2706*u3*u3+(774*u1-3687*u2)*u3+1356*u2*u2-579*u1*u2+64*u1*u1)/180.0;

//   if(pars.iStep<100) { beta0=0; beta1=0; beta2=0; }
  ftype2 w0 = d0*make_ftype2(1,1);
  ftype2 w1 = d1*make_ftype2(1,1);
  ftype2 w2 = d2*make_ftype2(1,1);
  
  ftype smooth_coeff=1e-6;
  w0.x/= (smooth_coeff+beta0.x)*(smooth_coeff+beta0.x); w0.y/= (smooth_coeff+beta0.y)*(smooth_coeff+beta0.y);
  w1.x/= (smooth_coeff+beta1.x)*(smooth_coeff+beta1.x); w1.y/= (smooth_coeff+beta1.y)*(smooth_coeff+beta1.y);
  w2.x/= (smooth_coeff+beta2.x)*(smooth_coeff+beta2.x); w2.y/= (smooth_coeff+beta2.y)*(smooth_coeff+beta2.y);
  
  ftype2 wsum = w0+w1+w2;
  w0/= wsum; w1/= wsum; w2/= wsum;

  //if(int(crd.y)==0 && int(crd.z==0) && fabs(crd.x+dir-1024)<=4 && i==0)
  //printf("i=%d niter=%d crd.x=%.10f, w012=%g %g %g, d012=%g,%g,%g, beta012=%g,%g,%g u123456=%g,%g,%g,%g,%g,%g S012=%g %g %g\n",i,niter,crd.x, w0,w1,w2, d0,d1,d2, beta0,beta1,beta2, u1,u2,u3,u4,u5,u6, S0,S1,S2);
  //printf("i=%d niter=%d crd.x=%.10f, w012=%g %g %g, S012=%g %g %g\n",i,niter,crd.x, w0,w1,w2, S0,S1,S2);

  return w0*S0+w1*S1+w2*S2;
}
inline __device__ ftype2 propagateNR_1D( ftype3 crd, int3 crd_cnt, const int i, const ftype4 base_gauge){
  #ifdef NRMESH
  const int3 crd_left = crd_cnt-make_int3(1,0,0);
  ftype2 fpr=make_ftype2(0,0);
  const int N=3;
  ftype bp[N];
  Cell c[N];
  for(int ix=0; ix<N; ix++) {
    int3 pos = crd_left+make_int3(ix,0,0);
    pos.x = (pos.x+pars.Nx)%pars.Nx;
    pos.y = (pos.y+pars.Ny)%pars.Ny;
    pos.z = (pos.z+pars.Nz)%pars.Nz;
    c[ix].load(pos.x,pos.y,pos.z);
    bp[ix] = c[ix].p.x-c[0].p.x;
    //kostyli
    if(bp[ix]<-pars.Nx/2) bp[ix]+=pars.Nx; 
  }
  ftype3 shift = crd-c[0].p;
  if(shift.x<-pars.Nx/2) shift.x+=pars.Nx;
  for(int ix=0; ix<N; ix++) {
    fpr+= LagrPolNR_1D<N>(ix,shift.x, bp)*c[ix].transfer_gauge(base_gauge, i);
    //int3 pos = crd_left+make_int3(ix,0,0);
    //printf("ix=%d pos=%d shift=%g coeff=%g\n",ix,pos.x,crd.x-c[0].p.x,LagrPolNR_1D<N>(ix,crd.x-bp[0], bp));
  }
  return fpr;
  #endif// NRMESH
}
inline __device__ ftype2 propagate3( ftype3 crd, int3 crd_left, const int i, const ftype4 base_gauge){
  const ftype3 shift = crd-make_ftype3(crd_left);
  ftype2 fpr=make_ftype2(0,0);
  const int N=3;
  for(int ix=0; ix<N; ix++) for(int iy=0; iy<N; iy++) for(int iz=0; iz<N; iz++) {
    int3 pos = crd_left+make_int3(ix,iy,iz);
    pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx;
    pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny;
    pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
    Cell c; c.load(pos.x,pos.y,pos.z);
    if(N==2)     fpr+= LinPol    (ix,iy,iz,shift.x,shift.y,shift.z)*c.transfer_gauge(base_gauge, i);
    else if(N>2) fpr+= LagrPol<N>(ix,iy,iz,shift.x,shift.y,shift.z)*c.transfer_gauge(base_gauge, i);
  }
  return fpr;
}
template<int N> inline __device__ ftype2 LagrEval( ftype3 crd, int3 crd_left, const int i, const ftype4 base_gauge){
  const ftype3 shift = crd-make_ftype3(crd_left);
  ftype2 fpr=make_ftype2(0,0);
  for(int ix=0; ix<N; ix++) for(int iy=0; iy<N; iy++) for(int iz=0; iz<N; iz++) {
    int3 pos = crd_left+make_int3(ix,iy,iz);
    pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx;
    pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny;
    pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
    Cell c; c.load(pos.x,pos.y,pos.z);
    if(N==2)     fpr+= LinPol    (ix,iy,iz,shift.x,shift.y,shift.z)*c.transfer_gauge(base_gauge, i);
    else if(N>2) fpr+= LagrPol<N>(ix,iy,iz,shift.x,shift.y,shift.z)*c.transfer_gauge(base_gauge, i);
  }
  return fpr;
}
template<int N> inline __device__ ftype3 LagrEvalMacro( ftype3 crd, int3 crd_left, const int i, const ftype4 base_gauge){
  const ftype3 shift = crd-make_ftype3(crd_left);
  ftype rho=0,mom=0,enrg=0;
  for(int ix=0; ix<N; ix++) for(int iy=0; iy<N; iy++) for(int iz=0; iz<N; iz++) {
    int3 pos = crd_left+make_int3(ix,iy,iz);
    pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx;
    pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny;
    pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
    Cell c; c.load(pos.x,pos.y,pos.z);
    ftype coeff = LagrPol<N>(ix,iy,iz,shift.x,shift.y,shift.z);
    rho += coeff*c.rho;
    mom += coeff*c.rho*c.uT.x;
    enrg+= coeff*c.rho*(Dim/2.0*c.uT.w+0.5*(c.uT.x*c.uT.x+c.uT.y*c.uT.y+c.uT.z*c.uT.z));
  }
  return make_ftype3(rho,mom,enrg);
  
}
inline __device__ ftype2 propagateTVD2( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir ){
  int3 crd_left;
  const int N=3;
  if(N==3)   crd_left = make_int3(floor(crd.x+0.5)-1, floor(crd.y+0.5), floor(crd.z+0.5));
  const ftype3 shift = crd-make_ftype3(crd_left);
  ftype2 fpr=make_ftype2(0,0);
  int part=0; //left part
  if(floor(crd.x+0.5)==floor(crd.x)) part=1; //right part
  int3 pos = crd_left+make_int3(0,0,0);
  pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx;
  pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny;
  pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
  Cell c; 
  ftype2 val0,val1,val2;
  ftype rhom,rhoc,rhop;
  c.load((pos.x+0)%pars.Nx,pos.y,pos.z); val0 = c.transfer_gauge(base_gauge, i, 0, dir); rhom=c.rho;
  c.load((pos.x+1)%pars.Nx,pos.y,pos.z); val1 = c.transfer_gauge(base_gauge, i, 1, dir); rhoc=c.rho;
  c.load((pos.x+2)%pars.Nx,pos.y,pos.z); val2 = c.transfer_gauge(base_gauge, i, 2, dir); rhop=c.rho;
  ftype coeff0 = (shift.x-1)*(shift.x-2)/2;
  ftype coeff1 = shift.x*(2-shift.x);
  ftype coeff2 = shift.x*(shift.x-1)/2;
  fpr = coeff0*val0 + coeff1*val1 + coeff2*val2;
  pars.cinfo[Cinfo::ind_zip(pos.x,pos.y,pos.z)].niter=0;
  if((rhom-rhoc)*(rhoc-rhop)<0) {
    pars.cinfo[Cinfo::ind_zip(pos.x,pos.y,pos.z)].niter=10;
    fpr = val1+(val2-val0)/2*(shift.x-1)+ 0.5*(val0+val2-2*val1)*(shift.x-1)*(shift.x-1);
    //if(part==0 && fabs(val1.x-val0.x)<fabs(val2.x-val1.x)) fpr = val1+(val1-val0)*(shift.x-1)+ 10*0.5*(val0+val2-2*val1)*(shift.x-1)*(shift.x-1);
    //if(part==1 && fabs(val1.x-val0.x)>fabs(val2.x-val1.x)) fpr = val1+(val2-val1)*(shift.x-1)+ 10*0.5*(val0+val2-2*val1)*(shift.x-1)*(shift.x-1);
  }
  return fpr;
}
template<int N> inline __device__ ftype2 propagate( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir ){
  int3 crd_left;
  //if(N==3 && dir<0)   crd_left = make_int3(floor(crd.x)-1, floor(crd.y+0.5)-1, floor(crd.z+0.5)-1); else
  //if(N==3 && dir>=0)   crd_left = make_int3(floor(crd.x), floor(crd.y+0.5)-1, floor(crd.z+0.5)-1);
  if(N==1 && dir<0)   crd_left = make_int3(floor(crd.x), floor(crd.y+0.5), floor(crd.z+0.5));
  if(N==1 && dir>=0)  crd_left = make_int3(floor(crd.x)+1, floor(crd.y+0.5), floor(crd.z+0.5));
  if(N==2)   crd_left = make_int3(floor(crd.x), floor(crd.y), floor(crd.z));
  if(N==3)   crd_left = make_int3(floor(crd.x+0.5)-1, floor(crd.y+0.5)-1, floor(crd.z+0.5)-1);
  if(N==4)   crd_left = make_int3(floor(crd.x)-1, floor(crd.y)-1, floor(crd.z)-1);
  if(N==5)   crd_left = make_int3(floor(crd.x+0.5)-2, floor(crd.y+0.5)-2, floor(crd.z+0.5)-2);
  if(N%2==0) crd_left = make_int3(floor(crd.x)-N/2+1, floor(crd.y)-N/2+1, floor(crd.z)-N/2+1);
  if(N%2==1) crd_left = make_int3(floor(crd.x+0.5)-(N-1)/2, floor(crd.y)-(N-1)/2, floor(crd.z)-(N-1)/2);
  const ftype3 shift = crd-make_ftype3(crd_left);
  ftype2 fpr=make_ftype2(0,0);

  for(int ix=0; ix<N; ix++) for(int iy=0; iy<N; iy++) for(int iz=0; iz<N; iz++) {
    int3 pos = crd_left+make_int3(ix,iy,iz);

    int doReverseX=0,doReverseY=0,doReverseZ=0;
    const bool outflow=0;
    if(outflow) {
      if(pos.x<0       ) { pos.x = -1-pos.x;                doReverseX=1; }
      if(pos.x>=pars.Nx) { pos.x = pars.Nx-1-pos.x+pars.Nx; doReverseX=1; }
      if(pos.y<0       ) { pos.y = -1-pos.y;                doReverseY=1; }
      if(pos.y>=pars.Ny) { pos.y = pars.Ny-1-pos.y+pars.Ny; doReverseY=1; }
      if(pos.z<0       ) { pos.z = -1-pos.z;                doReverseZ=1; }
      if(pos.z>=pars.Nz) { pos.z = pars.Nz-1-pos.z+pars.Nz; doReverseZ=1; }
    }
    
    pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx;
    pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny;
    pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
    Cell c; c.load(pos.x,pos.y,pos.z);

    if(doReverseX) c.mirrX();
    if(doReverseY) c.mirrY();
    if(doReverseZ) c.mirrZ();

    if(N==2)     fpr+= LinPol    (ix,iy,iz,shift.x,shift.y,shift.z)*c.transfer_gauge(base_gauge, i, ix, dir);
    else if(N>2) fpr+= LagrPol<N>(ix,iy,iz,shift.x,shift.y,shift.z)*c.transfer_gauge(base_gauge, i, ix, dir);
  }
  if(N==1){
    int3 pos = crd_left;
    pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx;
    pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny;
    pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
    Cell c; c.load(pos.x,pos.y,pos.z);
    return c.transfer_gauge(base_gauge, i);
  }
  return fpr;
}
template<int N> inline __device__ ftype2 propagateMacro( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir ){
  int3 crd_left;
  //if(N==3 && dir<0)   crd_left = make_int3(floor(crd.x)-1, floor(crd.y+0.5)-1, floor(crd.z+0.5)-1); else
  //if(N==3 && dir>=0)   crd_left = make_int3(floor(crd.x), floor(crd.y+0.5)-1, floor(crd.z+0.5)-1);
  if(N==2)   crd_left = make_int3(floor(crd.x), floor(crd.y), floor(crd.z));
  if(N==3)   crd_left = make_int3(floor(crd.x+0.5)-1, floor(crd.y+0.5)-1, floor(crd.z+0.5)-1);
  if(N==4)   crd_left = make_int3(floor(crd.x)-1, floor(crd.y)-1, floor(crd.z)-1);
  if(N==5)   crd_left = make_int3(floor(crd.x+0.5)-2, floor(crd.y+0.5)-2, floor(crd.z+0.5)-2);
  if(N%2==0) crd_left = make_int3(floor(crd.x)-N/2+1, floor(crd.y)-N/2+1, floor(crd.z)-N/2+1);
  const ftype3 shift = crd-make_ftype3(crd_left);
  ftype rho=0, mom=0, enrg=0;
  for(int ix=0; ix<N; ix++) for(int iy=0; iy<N; iy++) for(int iz=0; iz<N; iz++) {
    int3 pos = crd_left+make_int3(ix,iy,iz);
    pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx;
    pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny;
    pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
    Cell c; c.load(pos.x,pos.y,pos.z);
    ftype coeff = LagrPol<N>(ix,iy,iz,shift.x,shift.y,shift.z);
    rho += coeff*c.rho;
    mom += coeff*c.rho*c.uT.x;
    enrg+= coeff*c.rho*(2./Dim*c.uT.w+0.5*(c.uT.x*c.uT.x+c.uT.y*c.uT.y+c.uT.z*c.uT.z));
//     if(N==2)     fpr+= LinPol    (ix,iy,iz,shift.x,shift.y,shift.z)*c.transfer_gauge(base_gauge, i);
//     else if(N>2) fpr+= LagrPol<N>(ix,iy,iz,shift.x,shift.y,shift.z)*c.transfer_gauge(base_gauge, i);
  }
  ftype4 uT=make_ftype4(0,0,0,0);
  if(rho!=0) { uT.x= mom/rho; uT.w=enrg/rho-0.5*(uT.x*uT.x+uT.y*uT.y+uT.z*uT.z); uT.w*=2.0/Dim; }
  Cell cm; cm.set_base(rho, uT);
  return cm.transfer_gauge(base_gauge,i);
}
template<int N> inline __device__ ftype2 propagateWENOmacro( ftype3 crd, const int i, const ftype4 base_gauge, const ftype dir, const int niter){
  int3 crd_left;
  ftype3 S0 = LagrEvalMacro<4>(crd, make_int3(floor(crd.x)  , floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype3 S1 = LagrEvalMacro<4>(crd, make_int3(floor(crd.x)-1, floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype3 S2 = LagrEvalMacro<4>(crd, make_int3(floor(crd.x)-2, floor(crd.y+0.5), floor(crd.z+0.5)), i, base_gauge);
  ftype cdist = crd.x-(int(floor(crd.x)));
  ftype d0=(cdist*cdist+3*cdist+2)/20.;
  ftype d2=(cdist*cdist-5*cdist+6)/20.;
  ftype d1=1.-d0-d2;
  crd_left = make_int3(floor(crd.x)-2,floor(crd.x+0.5),floor(crd.x+0.5));
  int3 pos = crd_left+make_int3(0,0,0);
  int3 pos1 = crd_left+make_int3(0,0,0);
  int3 pos2 = crd_left+make_int3(1,0,0);
  int3 pos3 = crd_left+make_int3(2,0,0);
  int3 pos4 = crd_left+make_int3(3,0,0);
  int3 pos5 = crd_left+make_int3(4,0,0);
  int3 pos6 = crd_left+make_int3(5,0,0);
  pos.x = (pos.x%pars.Nx+pars.Nx)%pars.Nx; pos.y = (pos.y%pars.Ny+pars.Ny)%pars.Ny; pos.z = (pos.z%pars.Nz+pars.Nz)%pars.Nz;
  pos1.x= (pos1.x%pars.Nx+pars.Nx)%pars.Nx; pos1.y= (pos1.y%pars.Ny+pars.Ny)%pars.Ny; pos1.z = (pos1.z%pars.Nz+pars.Nz)%pars.Nz;
  pos2.x= (pos2.x%pars.Nx+pars.Nx)%pars.Nx; pos2.y= (pos2.y%pars.Ny+pars.Ny)%pars.Ny; pos2.z = (pos2.z%pars.Nz+pars.Nz)%pars.Nz;
  pos3.x= (pos3.x%pars.Nx+pars.Nx)%pars.Nx; pos3.y= (pos3.y%pars.Ny+pars.Ny)%pars.Ny; pos3.z = (pos3.z%pars.Nz+pars.Nz)%pars.Nz;
  pos4.x= (pos4.x%pars.Nx+pars.Nx)%pars.Nx; pos4.y= (pos4.y%pars.Ny+pars.Ny)%pars.Ny; pos4.z = (pos4.z%pars.Nz+pars.Nz)%pars.Nz;
  pos5.x= (pos5.x%pars.Nx+pars.Nx)%pars.Nx; pos5.y= (pos5.y%pars.Ny+pars.Ny)%pars.Ny; pos5.z = (pos5.z%pars.Nz+pars.Nz)%pars.Nz;
  pos6.x= (pos6.x%pars.Nx+pars.Nx)%pars.Nx; pos6.y= (pos6.y%pars.Ny+pars.Ny)%pars.Ny; pos6.z = (pos6.z%pars.Nz+pars.Nz)%pars.Nz;
  Cell c;
  c.load(pos1.x,pos.y,pos.z); ftype3 vals1 = make_ftype3(c.rho, c.rho*c.uT.x, c.rho*(Dim*c.uT.w+c.uT.x*c.uT.x+c.uT.y*c.uT.y+c.uT.z*c.uT.z));
  c.load(pos2.x,pos.y,pos.z); ftype3 vals2 = make_ftype3(c.rho, c.rho*c.uT.x, c.rho*(Dim*c.uT.w+c.uT.x*c.uT.x+c.uT.y*c.uT.y+c.uT.z*c.uT.z));
  c.load(pos3.x,pos.y,pos.z); ftype3 vals3 = make_ftype3(c.rho, c.rho*c.uT.x, c.rho*(Dim*c.uT.w+c.uT.x*c.uT.x+c.uT.y*c.uT.y+c.uT.z*c.uT.z));
  c.load(pos4.x,pos.y,pos.z); ftype3 vals4 = make_ftype3(c.rho, c.rho*c.uT.x, c.rho*(Dim*c.uT.w+c.uT.x*c.uT.x+c.uT.y*c.uT.y+c.uT.z*c.uT.z));
  c.load(pos5.x,pos.y,pos.z); ftype3 vals5 = make_ftype3(c.rho, c.rho*c.uT.x, c.rho*(Dim*c.uT.w+c.uT.x*c.uT.x+c.uT.y*c.uT.y+c.uT.z*c.uT.z));
  c.load(pos6.x,pos.y,pos.z); ftype3 vals6 = make_ftype3(c.rho, c.rho*c.uT.x, c.rho*(Dim*c.uT.w+c.uT.x*c.uT.x+c.uT.y*c.uT.y+c.uT.z*c.uT.z));

  ftype3 beta0 = (64*vals6*vals6+(-579*vals5+774*vals4-323*vals3)*vals6+1356*vals5*vals5+(1554*vals3-3687*vals4)*vals5+2706*vals4*vals4-2499*vals3*vals4+634*vals3*vals3)/180.0;
  ftype3 beta1 = (64*vals5*vals5+(-189*vals4-6*vals3+67*vals2)*vals5+366*vals4*vals4+(-537*vals3-6*vals2)*vals4+366*vals3*vals3-189*vals2*vals3+64*vals2*vals2)/180.0;
  ftype3 beta2 = (634*vals4*vals4+(-2499*vals3+1554*vals2-323*vals1)*vals4+2706*vals3*vals3+(774*vals1-3687*vals2)*vals3+1356*vals2*vals2-579*vals1*vals2+64*vals1*vals1)/180.0;

//   if(pars.iStep<100) { beta0=0; beta1=0; beta2=0; }
  ftype3 w0 = d0*make_ftype3(1,1,1);
  ftype3 w1 = d1*make_ftype3(1,1,1);
  ftype3 w2 = d2*make_ftype3(1,1,1);
  
  w0.x/= (1e-6+beta0.x)*(1e-6+beta0.x); w0.y/= (1e-6+beta0.y)*(1e-6+beta0.y); w0.z/= (1e-6+beta0.z)*(1e-6+beta0.z);
  w1.x/= (1e-6+beta1.x)*(1e-6+beta1.x); w1.y/= (1e-6+beta1.y)*(1e-6+beta1.y); w1.z/= (1e-6+beta1.z)*(1e-6+beta1.z);
  w2.x/= (1e-6+beta2.x)*(1e-6+beta2.x); w2.y/= (1e-6+beta2.y)*(1e-6+beta2.y); w2.z/= (1e-6+beta2.z)*(1e-6+beta2.z);
  
  ftype3 wsum = w0+w1+w2;
  w0/= wsum; w1/= wsum; w2/= wsum;
  
  ftype3 rme = w0*S0+w1*S1+w2*S2;
  ftype rho=rme.x, mom=rme.y, enrg=rme.z;
  
  ftype4 uT=make_ftype4(0,0,0,0);
  if(rho!=0) { uT.x= mom/rho; uT.w=enrg/rho-0.5*(uT.x*uT.x+uT.y*uT.y+uT.z*uT.z); uT.w*=2.0/Dim; }
  Cell cm; cm.set_base(rho, uT);
  return cm.transfer_gauge(base_gauge,i);
}
template<int N> inline __device__ ftype2 propagate_nonequal_1D( ftype3 crd, const int i, const ftype4 base_gauge, int niter ){
  int3 crd_left,crd_rght;
  if(N==2) {
    crd_left = make_int3(floor(crd.x), floor(crd.y), floor(crd.z));
    crd_rght = make_int3(crd_left.x+1, crd_left.y, crd_left.z);
    int3 posL = crd_left, posR = make_int3(posL.x+1, posL.y, posL.z);
    posL.x = (posL.x%pars.Nx+pars.Nx)%pars.Nx;    posR.x = (posR.x%pars.Nx+pars.Nx)%pars.Nx;
    posL.y = (posL.y%pars.Ny+pars.Ny)%pars.Ny;    posR.y = (posR.y%pars.Ny+pars.Ny)%pars.Ny;
    posL.z = (posL.z%pars.Nz+pars.Nz)%pars.Nz;    posR.z = (posR.z%pars.Nz+pars.Nz)%pars.Nz;
    ftype3 shiftL = crd-make_ftype3(crd_left), shiftR = make_ftype3(crd_rght)-crd;
    int cset = pars.cinfo[Cinfo::ind_zip(posL.x,posL.y,posL.z)].set;
    while(cset==0) {
      posL.x = (posL.x-1+pars.Nx)%pars.Nx; shiftL.x++;
      cset = pars.cinfo[Cinfo::ind_zip(posL.x,posL.y,posL.z)].set;
    }
    cset = pars.cinfo[Cinfo::ind_zip(posR.x,posR.y,posR.z)].set;
    while(cset==0) {
      posR.x = (posR.x+1+pars.Nx)%pars.Nx; shiftR.x++;
      cset = pars.cinfo[Cinfo::ind_zip(posR.x,posR.y,posR.z)].set;
    }
    ftype2 fpr=make_ftype2(0,0);
    Cell cl; cl.load(posL.x,posL.y,posL.z);
    Cell cr; cr.load(posR.x,posR.y,posR.z);
    fpr = shiftR.x/(shiftL.x+shiftR.x)*cl.transfer_gauge(base_gauge, i) +
          shiftL.x/(shiftL.x+shiftR.x)*cr.transfer_gauge(base_gauge, i);
    return fpr;
  }
}
template __device__ void LBMpull<2>(const int3 coord);
template __device__ void LBMpull<4>(const int3 coord);
template __device__ void LBMpull<6>(const int3 coord);
template __device__ void LBMpull<8>(const int3 coord);
