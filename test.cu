__device__ void test(float a, float b, float* c){
 
 //c[threadIdx.x]=a/b;
 c[threadIdx.x]=__fdiv_rd(a,b);
 //c[threadIdx.x]= __fdividef(a,b);
 //c[threadIdx.x]=(a+1)/2;

}
