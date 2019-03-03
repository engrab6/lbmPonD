#include "structs.cuh"

void Cell::set(const ftype r, const ftype4 vel){
  rho=r;
  uT=vel;
  ftype3 u3=make_ftype3(uT.x,uT.y,uT.z);
  const ftype mxw0 = ftype(1) - dot(u3,u3)*ftype(0.5)*dcs2;
  for(int i=0; i<Qn; i++) {
    ftype3 eidx = make_ftype3(e_c[i].x, e_c[i].y, e_c[i].z);
    ftype eu =  dot(eidx,u3)*dcs2;
    ftype mxw   = mxw0 + eu + eu*eu*ftype(0.5);
    f[i] = w_c[i]*rho*mxw;
  }
}
void Cell::set_base(const ftype r, const ftype4 vel){
  rho=r; uT=vel;
  for(int i=0; i<Qn; i++) f[i] = w_c[i]*rho;
  #ifdef DDF
  for(int i=0; i<Qn; i++) h[i] = w_c[i]*rho*(Dim*uT.w+uT.x*uT.x+uT.y*uT.y+uT.z*uT.z);
  #endif
}
void CellsSoA::set(const long int index, const ftype r, const ftype4 vel){
  rho[index]=r;
  ftype4 u=vel;
  ftype3 u3=make_ftype3(u.x,u.y,u.z);
  const ftype mxw0 = ftype(1) - dot(u3,u3)*ftype(0.5)*dcs2;
  for(int i=0; i<Qn; i++) {
    ftype3 eidx = make_ftype3(e_c[i].x, e_c[i].y, e_c[i].z);
    ftype eu =  dot(eidx,u3)*dcs2;
    ftype mxw   = mxw0 + eu + eu*eu*ftype(0.5);
    f[i][index] = w_c[i]*r*mxw;
  }
}
