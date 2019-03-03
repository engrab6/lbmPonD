#include "structs.cuh"
void drop(std::string* fprefix) {
  Group* gh = new Group[parsHost.Nx*parsHost.Ny*parsHost.Nz/Gsize/Gsize/Gsize];
  CHECK_ERROR(cudaMemcpy(gh, parsHost.loadGroups, parsHost.Nx*parsHost.Ny*parsHost.Nz/Gsize/Gsize/Gsize*sizeof(Group), cudaMemcpyDeviceToHost) );
 
  char fname1[256];
  FILE* pFile[5]; 
  sprintf(fname1, "%s/Rho_%07d.arr", fprefix->c_str(), parsHost.iStep); pFile[0] = fopen(fname1,"w");
  sprintf(fname1, "%s/Vx_%07d.arr" , fprefix->c_str(), parsHost.iStep); pFile[1] = fopen(fname1,"w");
  sprintf(fname1, "%s/Vy_%07d.arr" , fprefix->c_str(), parsHost.iStep); pFile[2] = fopen(fname1,"w");
  sprintf(fname1, "%s/Vz_%07d.arr" , fprefix->c_str(), parsHost.iStep); pFile[3] = fopen(fname1,"w");
  sprintf(fname1, "%s/T_%07d.arr"  , fprefix->c_str(), parsHost.iStep); pFile[4] = fopen(fname1,"w");
  int zero=0, twelve = 0*12, dim = 3, sizeofT = sizeof(ftype);  
  for(int i=0; i<sizeof(pFile)/sizeof(pFile[0]); i++) fwrite(&twelve     , sizeof(int), 1, pFile[i]);  //size of comment
  for(int i=0; i<sizeof(pFile)/sizeof(pFile[0]); i++) for(int ii=0; ii<twelve/4; ii++) fwrite(&zero   , sizeof(int), 1, pFile[i]);    // comment
  for(int i=0; i<sizeof(pFile)/sizeof(pFile[0]); i++) fwrite(&dim        , sizeof(int), 1, pFile[i]); //dim = 
  for(int i=0; i<sizeof(pFile)/sizeof(pFile[0]); i++) fwrite(&sizeofT    , sizeof(int), 1, pFile[i]); //data size
  for(int i=0; i<sizeof(pFile)/sizeof(pFile[0]); i++) fwrite(&parsHost.Nx, sizeof(int), 1, pFile[i]);
  for(int i=0; i<sizeof(pFile)/sizeof(pFile[0]); i++) fwrite(&parsHost.Ny, sizeof(int), 1, pFile[i]);
  for(int i=0; i<sizeof(pFile)/sizeof(pFile[0]); i++) fwrite(&parsHost.Nz, sizeof(int), 1, pFile[i]);
  //printf("saving %s\n",fname);
  for(int z=0; z<parsHost.Nz; z++) for(int y=0; y<parsHost.Ny; y++) for(int x=0; x<parsHost.Nx; x++) {
    const int2 gcind = Group::ind_conv(x, y, z);
    const int gindex=gcind.x, cindex=gcind.y;
    Cell c; gh[gindex].unpack(c, cindex);
    fwrite(&c.rho , sizeof(ftype), 1, pFile[0]);
    fwrite(&c.uT.x, sizeof(ftype), 1, pFile[1]);
    fwrite(&c.uT.y, sizeof(ftype), 1, pFile[2]);
    fwrite(&c.uT.z, sizeof(ftype), 1, pFile[3]);
    fwrite(&c.uT.w, sizeof(ftype), 1, pFile[4]);
  }
  delete[] gh;
  for(int i=0; i<sizeof(pFile)/sizeof(pFile[0]); i++) fclose(pFile[i]);
}
