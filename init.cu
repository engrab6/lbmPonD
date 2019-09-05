#include "structs.cuh"
__global__ void fill(){
  // Calculate surface coordinates
  unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
  unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
  unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;
  if (x < pars.Nx && y < pars.Ny && z<pars.Nz) {
    const ftype gamma = (Dim+2.0)/Dim;
    const ftype T0=1/10.;//*(gamma-1);//1./3./1000;//TLat;
    const ftype uinit = -1.0/10.;//*sqrt(T0);
    ftype rho0=3.0; ftype4 uT0=make_ftype4(uinit,0,0,T0);

    //if( (x-pars.Nx/2)*(x-pars.Nx/2)+(y-pars.Ny/2)*(y-pars.Ny/2)+(z-pars.Nz/2)*(z-pars.Nz/2)<10*10 ) {rho0=1.5; uT0=make_ftype4(0,0,0,TL);}
    if( (x>1024) ) {rho0=2.0; uT0=make_ftype4(-1.5/100.,0,0,0.25/10000./*0.75/10000.*/);}
    else {rho0=3.0; uT0=make_ftype4(-1.0/100.,0,0,/*1.0/10000.*/2./3./10000.);}
    if( (x>512) ) {rho0=1.0; uT0=make_ftype4(0.,0,0, 1.0/100000.);}
    else {rho0=1.0; uT0=make_ftype4(0.,0,0,100.0/100000.0);}
int xcnt=pars.Nx/2*0;
int ycnt=pars.Ny/2*0;
int zcnt=pars.Nz/2*0;
// ftype r = sqrt(ftype((x-xcnt)*(x-xcnt)+(y-ycnt)*(y-ycnt)+(z-zcnt)*(z-zcnt)));
ftype r=x;
int lc=0;
ftype Tleft=100000.0/1000000.;
ftype Trght=1.0     /1000000.;
    if( (x>=xcnt+lc/2) ) {rho0=1.0; uT0=make_ftype4(0,0,0, Trght);}
    else if(x>xcnt-lc/2){rho0=1.0; uT0=make_ftype4(0,0,0, Tleft+(x-xcnt+lc/2)*(Trght-Tleft)/lc);}
//    else if(x<100){rho0=1.0; uT0=make_ftype4(0.,0,0, (Trght+(Tleft-Trght)/100.*x));}
    else {rho0=1.0; uT0=make_ftype4(0,0,0,Tleft);}
//    if(x==xcnt) rho0=1.1; else rho0=1; uT0=make_ftype4(0,0,0,T0);

    if( (r<=pars.Nx/2) ) {rho0=1.0; uT0=make_ftype4(0,0,0, Tleft);}
    else {rho0=1.0; uT0=make_ftype4(0,0,0,Trght);}

    if( (x<=10) ) {rho0=2.2222222222222223; uT0=make_ftype4(1.375/10.,0,0, 1.996875/100.);}
//     if( (r<=pars.Nx/2) ) {rho0=1.3513513513513513; uT0=make_ftype4(0.6499999999999999/10.,0,0, 1.9425000000000001/100.);}
    else {rho0=1.0; uT0=make_ftype4(0,0,0,1/100.);}
    if((x-60)*(x-60)+(y-pars.Ny/2)*(y-pars.Ny/2)+(z-pars.Nz/2)*(z-pars.Nz/2)<40*40) {rho0=3.0; uT0=make_ftype4(0,0,0,1/3./100.);}

    ftype rb = sqrt(ftype((x-xcnt)*(x-xcnt)+(y-ycnt)*(y-ycnt)+(z-zcnt)*(z-zcnt)));
    if(rb<1) {rho0=1.0; uT0=make_ftype4(0,0,0,T0*100);} else {rho0=1.0; uT0=make_ftype4(0,0,0,T0);}

     //---Burger's wave-----//
     const ftype P0=1;
     const ftype Ampl=0.02;
     const ftype wlength = pars.Nx/1;
     ftype P=P0+Ampl*sin((x-xcnt)*2*M_PI/wlength); 
     ftype T = T0*pow(P,(gamma-1)/gamma);
     rho0=P/T;
     uT0=make_ftype4(Ampl*sin((x-xcnt)*2*M_PI/wlength)/P0*T0/sqrt(gamma*T0),0,0,T);
     
     //-----------------//
     
     //----3D-TGVortex----//
     /*const ftype M0=0.01;
     const ftype V0=M0*sqrt(gamma*T0);
     const ftype ux= V0*sin(x*2*M_PI/pars.Nx)*cos(y*2*M_PI/pars.Ny)*cos(z*2*M_PI/pars.Nz);
     const ftype uy=-V0*cos(x*2*M_PI/pars.Nx)*sin(y*2*M_PI/pars.Ny)*cos(z*2*M_PI/pars.Nz);
     const ftype uz=0;
     rho0=1.0;
     const ftype P0=rho0*T0;
     const ftype P = rho0*V0*V0*(1./(gamma*M0*M0) + 1./16*(cos(2*x*2*M_PI/pars.Nx)+cos(2*y*2*M_PI/pars.Ny))*(2+cos(2*z*2*M_PI/pars.Nz+2)));
     uT0 = make_ftype4(ux,uy,uz,P/rho0);
     ftype visc=(1/pars.pc.dtau-0.5)*T0;
     if(x==0 && y==0 && z==0) printf("Re=%g\n",V0/visc);*/
     //-----------------//
     
     //------Hot bubble------//
     /*ftype r2=(x-pars.Nx/2)*(x-pars.Nx/2)+(y-pars.Ny/2)*(y-pars.Ny/2)+(z-pars.Nz/2)*(z-pars.Nz/2);
     if(r2<1*1) {
       rho0=1.0/10.0;
       uT0=make_ftype4(0,0,0,T0*10);
     } else{ 
       rho0=1.0;
       uT0=make_ftype4(0,0,0,T0);
     }*/
     //----------------------//
    
    //-------------Sod problem----------//
//     if( x<pars.Nx/2 ) { rho0=1; uT0=make_ftype4(0,0,0,0.01/rho0);
//     } else {            rho0=0.125; uT0=make_ftype4(0,0,0,0.001/rho0);   }
//     if( x<pars.Nx/2 ) { rho0=1; uT0=make_ftype4(0,0,0,0.1/rho0); } 
//     else {              rho0=1; uT0=make_ftype4(0,0,0,0.02/rho0);   }
    //----------------------------------//
    
    //------Spherical Sod problem----------//
//     rb = sqrt(ftype(x*x+y*y+z*z));
//     if(rb<pars.Nx/2) {rho0=1.0; uT0=make_ftype4(0,0,0,0.05/rho0);} else {rho0=0.125; uT0=make_ftype4(0,0,0,0.005/rho0);}
    //-------------------------------------//

    //if( (x>128) ) {rho0=1.0012913189257875; uT0=make_ftype4((-63.218364928909956+30.0)/1000.,0,0, 1.7202785717223354/1000000.);}
    //else {rho0=2.0; uT0=make_ftype4((-31.5+30.0)/1000.,0,0,1000.0/1000000.0);}
    //if( (x>1024) ) {rho0=1.0; uT0=make_ftype4(+2./1000.,0,0,1.0/3000000.);}
    //else {rho0=1.0; uT0=make_ftype4(-2./1000.,0,0,1.0/3000000.);}
    //if( abs(x-pars.Nx/2)<1 ) {rho0=1; uT0=make_ftype4(uinit*(x-pars.Nx/2)/1.,0,0,T0*1);}
    Cell c; c.set_base(rho0, uT0);
    //Cell c; c.set_base(1.0, make_ftype4(0.1,0,0,T0));
    const int2 gcind = Group::ind_conv(x, y, z);
    const int gindex=gcind.x, cindex=gcind.y;
    #ifdef NRMESH
    c.p=make_ftype3(x+0.5-fabs(ftype(x)/pars.Nx-0.5),y,z);
    c.p=make_ftype3(x,y,z);
    #endif
    pars.loadGroups [gindex].pack(c, cindex);
    pars.storeGroups[gindex].pack(c, cindex);
    if(x%1==0) pars.cinfo[Cinfo::ind_zip(x,y,z)].set=1;
  }
}
void init(){
  dim3 dimBlock(16, 16, 1);
  dim3 dimGrid((parsHost.Nx+dimBlock.x-1)/dimBlock.x, (parsHost.Ny+dimBlock.y-1)/dimBlock.y, (parsHost.Nz+dimBlock.z-1)/dimBlock.z);
    
  cuTimer t0;
  printf("LBM memory %g GB\n",parsHost.Nx*parsHost.Ny*parsHost.Nz*sizeof(Cell)/1024./1024./1024.);
  
//  CHECK_ERROR( cudaMalloc((void**)&parsHost.cells, (long int)parsHost.Nx*parsHost.Ny*parsHost.Nz*sizeof(Cell)) );
  
//    for(int i=0; i<Qn; i++) CHECK_ERROR( cudaMalloc((void**)&parsHost.csoa.f[i], (long int)parsHost.Nx*parsHost.Ny*parsHost.Nz*sizeof(ftype)) );
//    CHECK_ERROR( cudaMalloc((void**)&parsHost.csoa.rho, (long int)parsHost.Nx*parsHost.Ny*parsHost.Nz*sizeof(ftype)) );
  
  CHECK_ERROR( cudaMalloc((void**)&parsHost.loadGroups , (long int)parsHost.Nx*parsHost.Ny*parsHost.Nz/Gsize/Gsize/Gsize*sizeof(Group)) );
  CHECK_ERROR( cudaMalloc((void**)&parsHost.storeGroups, (long int)parsHost.Nx*parsHost.Ny*parsHost.Nz/Gsize/Gsize/Gsize*sizeof(Group)) );
  CHECK_ERROR( cudaMalloc((void**)&parsHost.cinfo, (long int)parsHost.Nx*parsHost.Ny*parsHost.Nz*sizeof(Cinfo)) );
  CHECK_ERROR( cudaMemset(parsHost.cinfo, 0, (long int)parsHost.Nx*parsHost.Ny*parsHost.Nz*sizeof(Cinfo)) );

  copy2dev( parsHost, pars );
  fill<<<dimGrid, dimBlock>>>(); cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() ); 
}

void LBMParamsHost::reset() {
  CHECK_ERROR( cudaMemset(cinfo, 0, (long int)Nx*Ny*Nz*sizeof(Cinfo)) );
}
