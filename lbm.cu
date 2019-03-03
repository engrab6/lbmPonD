//#include <cuda_fp16.h>
#include "cuda_math.h"
//#include "cuda_math_double.h"
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "structs.cuh"

#include "im2D.h"
#include "im3D.hpp"
int type_diag_flag=0;

im3D_pars im3DHost;
LBMParamsHost parsHost;
__constant__ LBMParams pars;
__constant__ int3 e[Qn];
__constant__ ftype3 ef[Qn];
__constant__ ftype w[Qn];
__constant__ int reverseX[Qn];  
__constant__ int reverseY[Qn];  
__constant__ int reverseZ[Qn];  
__constant__ int reverseXYZ[Qn];

const char* FuncStr[] = {"Rho","Velocity","Mach","T","M4mom","Iterations","f0","f1","f2","f3","f4","f5","f6","f7","f8"}; 
//const char* FuncStr[] = {"F","G"}; 

__device__ float pow2(float v) { return v*v; }
__global__ void lbm_draw_surf(float* buf){
}
__global__ void lbm_draw(float* buf) {
  for(int ix=threadIdx.x; ix<pars.Nx; ix+=pars.Nx/4) {
    int iz=blockIdx.x;
    int iy=blockIdx.y;
    int xbuf=ix;
    int ybuf=iy;
    float* pbuf=&buf[xbuf+pars.Nx*(ybuf+pars.Ny*iz)];
    ftype f[Qn];
    Cell cell;
    const int2 gcind = Group::ind_conv(ix, iy, iz);
    const int gindex=gcind.x, cindex=gcind.y;
    pars.loadGroups[gindex].unpack(cell, cindex);
    for(int i=0; i<Qn; i++) f[i] = cell.f[i];
    ftype rho=cell.rho;
    ftype4 uT=cell.uT;
    uT/=pars.dt; uT.w/=pars.dt;
    ftype M4mom=0;
    for(int i=0; i<Qn; i++) {
      ftype viX = sqrt(uT.w/TLat)*ef[i].x+uT.x;
      M4mom+= f[i]*viX*viX*viX*viX;
    }
    Cinfo& cinf = pars.cinfo[Cinfo::ind_zip(ix,iy,iz)];
    int niter=cinf.niter;
    if(cinf.set==0) { *pbuf=0; continue; }
    switch(pars.nFunc) {
      case 0: *pbuf=cell.rho; break;
      //case 1: *pbuf=length(make_ftype3(cell.uT.x,cell.uT.y,cell.uT.z)); break; 
      case 1: *pbuf=uT.x; break;
      case 2: *pbuf=uT.x/sqrt(uT.w); break; 
      case 3: *pbuf=uT.w; break;
      case 4: *pbuf=M4mom; break;
      case 5: *pbuf=float(niter); break; 
      case 6: *pbuf=f[0]; break;
      case 7: *pbuf=f[1]; break;
      case 8: *pbuf=f[2]; break;
      case 9: *pbuf=f[3]; break;
      case 10: *pbuf=f[4]; break;
      case 11: *pbuf=f[5]; break;
      case 12: *pbuf=f[6]; break;
      case 13: *pbuf=f[7]; break;
      case 14: *pbuf=f[8]; break;
    }
  }
}
void draw_all(){
  lbm_draw<<<dim3(parsHost.Nz,parsHost.Ny),parsHost.Nx/4>>>(parsHost.arr4im.Arr3Dbuf);
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  im3DHost.initCuda(parsHost.arr4im);
}

void init();
void drop(std::string* fprefix);
void calcStep();

void idle_func_calc::step() {
  for(int i=0;i<1; i++) calcStep(); im3DHost.save_png(parsHost.iStep);
  draw_all();
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  recalc_at_once=true;
}

static void key_func(unsigned char key, int x, int y) {
  if(type_diag_flag>=2) printf("keyN=%d, coors=(%d,%d)\n", key, x, y);
  if(key == 'h') {
    printf("\
======= Управление mxw3D:\n\
  <¦>  \tИзменение функции для визуализации: WEH¦Sx¦Ez¦Ey¦Ex¦Hx¦Hy¦Hz¦Sy¦Sz¦eps\n\
«Enter»\tПересчёт одного большого шага\n\
   b   \tвключает пересчёт в динамике (см. «Управление динамикой»)\n\
"); im3DHost.print_help();
    return;
  }
  ftype t0;
  switch(key) {
  //case '>': if(parsHost.nFunc<parsHost.MaxFunc) parsHost.nFunc++; break;
  //case '<': if(parsHost.nFunc>0) parsHost.nFunc--; break;
  case '>': parsHost.nFunc = (parsHost.nFunc+1)%parsHost.MaxFunc; break;
  case '<': parsHost.nFunc = (parsHost.nFunc+parsHost.MaxFunc-1)%parsHost.MaxFunc; break;
  case 13: for(int i=0;i<1; i++) calcStep(); break;
  default: if(!im3DHost.key_func(key, x, y)) {
  if(type_diag_flag>=0) printf("По клавише %d в позиции (%d,%d) нет никакого действия\n", key, x, y);
  } return;
  }
  copy2dev( parsHost, pars );
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  draw_all();
  recalc_at_once=true;
}
static void draw_func() {
  im3DHost.fName = FuncStr[parsHost.nFunc]; 
  glutPostRedisplay();
  im2D.draw(im3DHost.reset_title()); 
}

//void (*idle_func_ptr)(float* );
static void idle_func() { im3DHost.recalc_func(); }
static void mouse_func(int button, int state, int x, int y) { im3DHost.mouse_func(button, state, x, y); }
static void motion_func(int x, int y) { im3DHost.motion_func(x, y); }
static void special_func(int key, int x, int y) { 
  im3DHost.special_func(key, x, y);
}

bool interactive=true, test_only=false;
//bool help_only=false, test_only=false;
void LBMParamsHost::set(){
  nFunc = 0; MaxFunc = sizeof(FuncStr)/sizeof(char*);
  Nx=(1<<Rank)*brick.x; //128;
  Ny=(1<<Rank)*brick.y; //128;
  Nz=(1<<Rank)*brick.z; //128;
  Nx=LSizeX*brick.x; Ny=LSizeY*brick.y; Nz=LSizeZ*brick.z;
  if(!test_only) {
  cudaChannelFormatDesc channelDesc = cudaCreateChannelDesc(32, 32, 0, 0, cudaChannelFormatKindFloat);
  CHECK_ERROR( cudaGetLastError() );
  CHECK_ERROR( cudaMalloc3DArray(&data, &channelDesc, make_cudaExtent(Nx*No,Ny,Nz), cudaArraySurfaceLoadStore) );
  struct cudaResourceDesc resDesc; memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType = cudaResourceTypeArray;
  resDesc.res.array.array = data;
  texSdata=0;
  CHECK_ERROR( cudaCreateSurfaceObject(&texSdata, &resDesc) );
  
  channelDesc = cudaCreateChannelDesc(32, 0, 0, 0, cudaChannelFormatKindFloat);
  CHECK_ERROR( cudaGetLastError() );
  CHECK_ERROR( cudaMalloc3DArray(&fdata, &channelDesc, make_cudaExtent(Nx,Ny,Nz), cudaArraySurfaceLoadStore) );
  CHECK_ERROR( cudaMalloc3DArray(&gdata, &channelDesc, make_cudaExtent(Nx,Ny,Nz), cudaArraySurfaceLoadStore) );
  memset(&resDesc, 0, sizeof(resDesc));
  resDesc.resType = cudaResourceTypeArray;
  texFdata=0; texGdata=0;
  resDesc.res.array.array = fdata; CHECK_ERROR( cudaCreateSurfaceObject(&texFdata, &resDesc) );
  resDesc.res.array.array = gdata; CHECK_ERROR( cudaCreateSurfaceObject(&texGdata, &resDesc) );
  }

  CHECK_ERROR(cudaMallocManaged((void**)&NiterMax,sizeof(int   ))); CHECK_ERROR(cudaMemset(NiterMax, 0, sizeof(int  )));
  CHECK_ERROR(cudaMallocManaged((void**)&mass    ,sizeof(ftype ))); CHECK_ERROR(cudaMemset(mass    , 0, sizeof(ftype)));
  CHECK_ERROR(cudaMallocManaged((void**)&enrg    ,sizeof(ftype ))); CHECK_ERROR(cudaMemset(enrg    , 0, sizeof(ftype)));
  CHECK_ERROR(cudaMallocManaged((void**)&moment  ,sizeof(ftype3))); CHECK_ERROR(cudaMemset(moment  , 0, sizeof(ftype3)));
  dt=1;
  
  pc.set();
}

int print_help() {
  printf("help | using in test|batch mode:\n ./lbm [--help|--test|--batch]\n");
  printf("using in interactive mode:\n ./lbm %s\n", im3DHost.command_line_help_string());
  im3DHost.print_command_line_help();
  return 0;
}

void reset(im3D_pars* p=0);
int maxR;
int _main(int argc, char** argv) {
  maxR=1;//atoi(argv[1]);
  ::reset();
  argv ++; argc --;
  im3DHost.reset();
  while(argc>0 && strncmp(*argv,"--",2)==0) {
    int pp=1;
    if(strcmp(*argv,"--test")==0) test_only = true;
    else if(strcmp(*argv,"--batch")==0) interactive = false;
    else pp = im3DHost.init_from_command_line(argv);
    if(pp<=0) return print_help();
    else if(pp==1) printf("par: %s; \n", argv[0]);
    else if(pp==2) printf("par: %s; vals: %s\n", argv[0], argv[1]);
    argv += pp; argc -= pp;
  };
  if(test_only) printf("No GL\n");
  else printf("With GL\n");
  im2D.get_device(3,0);
  type_diag_flag = 1;
try {
  if(type_diag_flag>=1) printf("Настройка опций визуализации по умолчанию\n");
  cudaTimer tm; tm.start();
  copy2dev( e_host[0], e ); copy2dev( w_host[0], w ); 
  #if defined D1Q3 || defined D3Q27
  #else
  copy2dev( ef_host[0], ef ); 
  #endif
  #ifdef D3Q125
  int reverseX_correct[Qn], reverseY_correct[Qn], reverseZ_correct[Qn];
  for(int i=0; i<Qn; i++) {
    reverseX_correct[i] = reverseX_host[i];
    reverseY_correct[i] = reverseY_host[i];
    reverseZ_correct[i] = reverseZ_host[i];
  }
  for(int i=0; i<Qn; i++) {
    if(ef_host[i].x!=0) {
      int j=0, found=0;
      while(found==0 && j<Qn) {
	if(ef_host[j]==ef_host[i]*make_ftype3(-1,1,1)) {found=1; break;}
	j++;
      }
      if(found==0) { printf("cannot find reverse index: it's strange\n"); exit(-1); }
      reverseX_correct[i] = j;
    }
    if(ef_host[i].y!=0) {
      int j=0, found=0;
      while(found==0 && j<Qn) {
	if(ef_host[j]==ef_host[i]*make_ftype3(1,-1,1)) {found=1; break;}
	j++;
      }
      if(found==0) { printf("cannot find reverse index: it's strange\n"); exit(-1); }
      reverseY_correct[i] = j;
    }
    if(ef_host[i].z!=0) {
      int j=0, found=0;
      while(found==0 && j<Qn) {
	if(ef_host[j]==ef_host[i]*make_ftype3(1,1,-1)) {found=1; break;}
	j++;
      }
      if(found==0) { printf("cannot find reverse index: it's strange\n"); exit(-1); }
      reverseZ_correct[i] = j;
    }
  }
  copy2dev( reverseX_correct[0], reverseX ); 
  copy2dev( reverseY_correct[0], reverseY ); 
  copy2dev( reverseZ_correct[0], reverseZ ); 
  #else
  copy2dev( reverseX_host[0], reverseX ); 
  copy2dev( reverseY_host[0], reverseY ); 
  copy2dev( reverseZ_host[0], reverseZ ); 
  #endif
  copy2dev( reverseXYZ_host[0], reverseXYZ ); 
  parsHost.set();
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  copy2dev( parsHost, pars );
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  init();
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );
  copy2dev( parsHost, pars );
  cudaDeviceSynchronize(); CHECK_ERROR( cudaGetLastError() );

  if(test_only) {
//     while(parsHost.iStep<1500*4*8) {
    while(parsHost.iStep<1) {
      const int Nx=parsHost.Nx, Ny=parsHost.Ny, Nz=parsHost.Nz;
      tm.start();
      calcStep();
      double tCpu=tm.stop();
      printf("run time: %.2f msec, %.2f MLU/sec\n", tCpu, 1.e-6*Nx*Ny*Nz/tCpu);
//       break;
    }
    return 0;
  }

  tm.start();
  parsHost.reset_im();
  im3DHost.reset(parsHost.arr4im);
  copy2dev( parsHost, pars );
  //CHECK_ERROR( cudaMemset(parsHost.arr4im.Arr3Dbuf,0,((long long int)Nx)*Ny*Nz*sizeof(ftype)) );
  im2D.get_device(3,0);
  im2D.init_image(argc,argv, im3DHost.bNx, im3DHost.bNy, "im3D");
  im3DHost.init3D(parsHost.arr4im); im3DHost.iz0=parsHost.Nx-1; im3DHost.key_func('b',0,0);
  im3DHost.initCuda(parsHost.arr4im);
  draw_all();

  if(type_diag_flag>=1) printf("Настройка GLUT и запуск интерфейса\n");
  glutIdleFunc(idle_func);
  glutKeyboardFunc(key_func);
  glutMouseFunc(mouse_func);
  glutMotionFunc(motion_func);
  glutDisplayFunc(draw_func);
  glutSpecialFunc(special_func);
  if(type_diag_flag>=0) printf("Init cuda device: %.1f msec\n", tm.stop());
  glutMainLoop();
} catch(...) {
  printf("Возникла какая-то ошибка.\n");
}
  parsHost.clear();
  return -1;
}
int main(int argc, char** argv) {
  return _main(argc,argv);
}
int run(int argc, char** argv) {
  return _main(argc,argv);
}
float get_val_from_arr3D(int ix, int iy, int iz) {
  Arr3D_pars& arr=parsHost.arr4im;
  if(arr.inCPUmem) return arr.Arr3Dbuf[arr.get_ind(ix,iy,iz)];
  float res=0.0;
  if(arr.inGPUmem) CHECK_ERROR(cudaMemcpy(&res, arr.get_ptr(ix,iy,iz), sizeof(float), cudaMemcpyDeviceToHost));
  return res;
}
