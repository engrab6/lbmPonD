__global__ __launch_bounds__(brick.x*brick.y*brick.z) void AMR_add(){
  const int3 coord = make_int3(blockIdx*brick+threadIdx);
  int cset = pars.cinfo[Cinfo::ind_zip(coord.x,coord.y,coord.z)].set;
  if(cset==0) return;
  Cell cm; cm.load_st(coord.x,coord.y,coord.z);
  int3 coordp=coord;
  int dist=0;
  do{
    coordp.x = (coordp.x+1)%pars.Nx;
    cset = pars.cinfo[Cinfo::ind_zip(coordp.x,coordp.y,coordp.z)].set;
    dist++;
  } while (cset==0);
  Cell cp; cp.load_st(coordp.x,coordp.y,coordp.z);
  ftype vdif = cp.uT.x-cm.uT.x;
  int3 coord_mid=coord; coord_mid.x=(coord_mid.x+dist/2)%pars.Nx;
  if(vdif*vdif/min(cm.uT.w,cp.uT.w)>0.01 || fabs(1-cm.rho/cp.rho)>0.01) {
    if(dist<=1) printf("limit of cells division exceeded\n");
    if(dist%2==0) {
      Cell c0;
      c0.rho = 0.5*(cm.rho+cp.rho);
      c0.uT = 0.5*(cm.uT+cp.uT);
      for(int i=0; i<Qn; i++) c0.f[i] = 0.5*(cm.f[i]+cp.f[i]);
      c0.save_12(coord_mid.x,coord_mid.y,coord_mid.z);
      pars.cinfo[Cinfo::ind_zip(coord_mid.x,coord_mid.y,coord_mid.z)].setnew=1;
    }
  }
}
__global__ __launch_bounds__(1) void AMR_remove(){
  for(int ix=0;ix<pars.Nx; ix++) {
    const int3 coord = make_int3(ix,0,0);
    int cset = pars.cinfo[Cinfo::ind_zip(coord.x,coord.y,coord.z)].set;
    if(cset==0) continue;
    Cell c0; c0.load_st(coord.x,coord.y,coord.z);
    int3 coordp=coord, coordm=coord;
    int distM=0,distP=0;
    do{
      coordm.x = (coordm.x-1+pars.Nx)%pars.Nx;
      cset = pars.cinfo[Cinfo::ind_zip(coordm.x,coordm.y,coordm.z)].set;
      distM++;
    } while (cset==0);
    do{
      coordp.x = (coordp.x+1)%pars.Nx;
      cset = pars.cinfo[Cinfo::ind_zip(coordp.x,coordp.y,coordp.z)].set;
      distP++;
    } while (cset==0);
    Cell cm; cm.load_st(coordm.x,coordm.y,coordm.z);
    Cell cp; cp.load_st(coordp.x,coordp.y,coordp.z);
    ftype vdifm = cm.uT.x-c0.uT.x;
    ftype vdifp = cp.uT.x-c0.uT.x;
    if(vdifm*vdifm/c0.uT.w<0.01 && vdifp*vdifp/c0.uT.w<0.01 && fabs(cp.rho/c0.rho-1)<0.01 && fabs(cm.rho/c0.rho-1)<0.01 && (distM<4 && distP<4)) {
      pars.cinfo[Cinfo::ind_zip(coord.x,coord.y,coord.z)].set=0;
      pars.cinfo[Cinfo::ind_zip(coord.x,coord.y,coord.z)].set=0;
    }
  }
}
