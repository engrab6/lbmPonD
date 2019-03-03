const ftype dx=1.0;
const ftype dy=1.0;
const ftype dz=1.0;
const ftype dt=1.0;

//#define DDF
  #define D1Q5
//   #define D1Q3
//#define D3Q125
// #define D3Q27

// #define NRMESH

#ifdef D3Q27
const int Qn=27;

const ftype TLat=1./3.;///100.;//1/3.;
const ftype Dim=3;

const int3 e_host[Qn] = {
 /*0 -0 */ { 0, 0, 0},
 /*1 -7 */ { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}, { 1, 1, 0}, { 1, 0, 1}, {0, 1, 1}, {1, 1, 1},
 /*8 -11*/ {-1, 0, 0}, {-1, 1, 0}, {-1, 0, 1}, {-1, 1, 1},
 /*12-15*/ { 0,-1, 0}, { 0,-1, 1}, { 1,-1, 0}, { 1,-1, 1}, 
 /*16-19*/ { 0, 0,-1}, { 1, 0,-1}, { 0, 1,-1}, { 1, 1,-1},
 /*20-21*/ {-1,-1, 0}, {-1,-1, 1},
 /*22-23*/ {-1, 0,-1}, {-1, 1,-1}, 
 /*24-25*/ { 0,-1,-1}, { 1,-1,-1},
 /*26-26*/ {-1,-1,-1}
};
#define E_MATRIX_X(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18) 0
#define E_MATRIX_Y(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18) 0
#define E_MATRIX_Z(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18) 0
const ftype W0=(1-TLat)*(1-TLat)*(1-TLat);
const ftype W1=(1-TLat)*(1-TLat)*TLat/2.0;
const ftype W2=(1-TLat)*TLat*TLat/4.0;
const ftype W3=TLat*TLat*TLat/8.;
const ftype w_host[Qn] = {W0,
  W1, W1, W1, W2, W2, W2, W3,
  W1, W2, W2, W3, 
  W1, W2, W2, W3, 
  W1, W2, W2, W3, 
  W2, W3,
  W2, W3,
  W2, W3,
  W3
};
/*const ftype W0=8./27;
const ftype W1=1./27;
const ftype W2=1./54;
const ftype W3=1./216;
const ftype w_host[Qn] = {8./27,
  2./27, 2./27, 2./27, 1./54, 1./54, 1./54, 1./216,
  2./27, 1./54, 1./54, 1./216, 
  2./27, 1./54, 1./54, 1./216, 
  2./27, 1./54, 1./54, 1./216, 
  1./54, 1./216,
  1./54, 1./216,
  1./54, 1./216,
  1./216
};*/
const int reverseX_host[Qn] = {0, 8, 2, 3, 9, 10, 6, 11, 1, 4, 5, 7, 12, 13, 20, 21, 16, 22, 18, 23, 14, 15, 17, 19, 24, 26, 25};
const int reverseY_host[Qn] = {0, 1, 12, 3, 14, 5, 13, 15, 8, 20, 10, 21, 2, 6, 4, 7, 16, 17, 24, 25, 9, 11, 22, 26, 18, 19, 23};
const int reverseZ_host[Qn] = {0, 1, 2, 16, 4, 17, 18, 19, 8, 9, 22, 23, 12, 24, 14, 25, 3, 5, 6, 7, 20, 26, 10, 11, 13, 15, 21};
const int reverseXYZ_host[Qn]={0, 8, 12, 16, 20, 22, 24, 26, 1, 14, 17, 25, 2, 18, 9, 23, 3, 10, 13, 21, 4, 19, 5, 15, 6, 11, 7};

#elif defined D3Q19
const int Qn=19;

const int3 e_host[Qn] = {
 { 0, 0, 0},
 { -1, 0, 0}, { 1, 0, 0}, 
 { -1,-1, 0}, { 0,-1, 0}, { 1,-1, 0}, {-1, 1, 0}, { 0, 1, 0}, { 1, 1, 0},
 { -1, 0,-1}, { 0, 0,-1}, { 1, 0,-1}, {-1, 0, 1}, { 0, 0, 1}, { 1, 0, 1}, 
 {  0,-1,-1},             { 0, 1,-1}, { 0,-1, 1},             { 0, 1, 1}
};
#define E_MATRIX_X(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18) -f1 +f2-f3 +f5-f6   +f8-f9+f11-f12+f14
#define E_MATRIX_Y(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18) -f3 -f4-f5+f6+f7+f8 -f15+f16-f17+f18
#define E_MATRIX_Z(f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16,f17,f18) -f9-f10-f11+f12+f13+f14-f15-f16+f17+f18
const ftype W0=1./3;
const ftype W1=1./18;
const ftype W2=1./36;
const ftype w_host[Qn] = {1./3,
  1./18, 1./18, 
  1./36, 1./18, 1./36, 1./36, 1./18, 1./36,
  1./36, 1./18, 1./36, 1./36, 1./18, 1./36, 
  1./36,        1./36, 1./36,        1./36
};
const int reverseX_host[Qn] = {0, 2, 1, 5, 4, 3, 8, 7, 6, 11, 10,  9, 14, 13, 12, 15, 16, 17, 18};
const int reverseY_host[Qn] = {0, 1, 2, 6, 7, 8, 3, 4, 5,  9, 10, 11, 12, 13, 14, 16, 15, 18, 17};
const int reverseZ_host[Qn] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14,  9, 10, 11, 17, 18, 15, 16};
const int reverseXYZ_host[Qn]={0, 2, 1, 8, 7, 6, 5, 4, 3, 14, 13, 12, 11, 10,  9, 18, 17, 16, 15};

#elif defined D1Q3
const int Qn=3;

const ftype TLat=1./3.;///100.;//1/3.;
const ftype Dim=1;

const int3 e_host[Qn] = {
 /*0 -0 */ { 0, 0, 0}, { 1, 0, 0}, { -1, 0, 0} 
};
const ftype W0=2./3.;
const ftype W1=1./6.;
const ftype w_host[Qn] = {W0,W1,W1};
const int reverseX_host[Qn] = {0, 2, 1};
const int reverseY_host[Qn] = {0, 1, 2};
const int reverseZ_host[Qn] = {0, 1, 2};
const int reverseXYZ_host[Qn]={0, 2, 1};

#elif defined D1Q5
const int Qn=5;

const ftype TLat=1;//(1.-0.6324555320336759);///100.;//1/3.;
const ftype Dim=1;

const ftype ec1=1.3556261799742657;
const ftype ec2=2.8569700138728056;
const ftype ec1d=ec1*ec1;
const ftype ec2d=ec2*ec2;
const int3 e_host[Qn] = {
 /*0 -0 */ { 0, 0, 0}, { int(ec1), 0, 0}, { -int(ec1), 0, 0}, {int(ec2), 0, 0}, {-int(ec2), 0, 0}
};
const ftype3 ef_host[Qn] = {
 /*0 -0 */ { 0, 0, 0}, { ec1, 0, 0}, { -ec1, 0, 0}, {ec2, 0, 0}, {-ec2, 0, 0}
};
const ftype W0=(3*TLat*TLat-(ec1d+ec2d)*TLat+ec1d*ec2d)/(ec1d*ec2d);
const ftype W1=TLat*(ec2d-3*TLat)/(2*ec1d*(ec2d-ec1d));
const ftype W2=TLat*(ec1d-3*TLat)/(2*ec2d*(ec1d-ec2d));
const ftype w_host[Qn] = {W0,W1,W1,W2,W2};
const int reverseX_host[Qn] = {0, 2, 1, 4,3};
const int reverseY_host[Qn] = {0, 1, 2,3,4};
const int reverseZ_host[Qn] = {0, 1, 2,3,4};
const int reverseXYZ_host[Qn]={0, 2, 1,4,3};

#elif defined D3Q125
const int Qn=125;

const ftype TLat=1;//(1.-0.6324555320336759);///100.;//1/3.;
const ftype Dim=3;

const ftype ec1=1.3556261799742657;
const ftype ec2=2.8569700138728056;
const ftype ec1d=ec1*ec1;
const ftype ec2d=ec2*ec2;
const int3 e_host[Qn] = {
  {0        , 0        , 0 }, {0        , 0        , +int(ec1) }, {0        , 0        , -int(ec1)}, {0        , 0        , +int(ec2)}, {0        , 0        , -int(ec2)},
  {0        , +int(ec1), 0 }, {0        , +int(ec1), +int(ec1) }, {0        , +int(ec1), -int(ec1)}, {0        , +int(ec1), +int(ec2)}, {0        , +int(ec1), -int(ec2)},
  {0        , -int(ec1), 0 }, {0        , -int(ec1), +int(ec1) }, {0        , -int(ec1), -int(ec1)}, {0        , -int(ec1), +int(ec2)}, {0        , -int(ec1), -int(ec2)},
  {0        , +int(ec2), 0 }, {0        , +int(ec2), +int(ec1) }, {0        , +int(ec2), -int(ec1)}, {0        , +int(ec2), +int(ec2)}, {0        , +int(ec2), -int(ec2)},
  {0        , -int(ec2), 0 }, {0        , -int(ec2), +int(ec1) }, {0        , -int(ec2), -int(ec1)}, {0        , -int(ec2), +int(ec2)}, {0        , -int(ec2), -int(ec2)},
  {+int(ec1), 0        , 0 }, {+int(ec1), 0        , +int(ec1) }, {+int(ec1), 0        , -int(ec1)}, {+int(ec1), 0        , +int(ec2)}, {+int(ec1), 0        , -int(ec2)},
  {+int(ec1), +int(ec1), 0 }, {+int(ec1), +int(ec1), +int(ec1) }, {+int(ec1), +int(ec1), -int(ec1)}, {+int(ec1), +int(ec1), +int(ec2)}, {+int(ec1), +int(ec1), -int(ec2)},
  {+int(ec1), -int(ec1), 0 }, {+int(ec1), -int(ec1), +int(ec1) }, {+int(ec1), -int(ec1), -int(ec1)}, {+int(ec1), -int(ec1), +int(ec2)}, {+int(ec1), -int(ec1), -int(ec2)},
  {+int(ec1), +int(ec2), 0 }, {+int(ec1), +int(ec2), +int(ec1) }, {+int(ec1), +int(ec2), -int(ec1)}, {+int(ec1), +int(ec2), +int(ec2)}, {+int(ec1), +int(ec2), -int(ec2)},
  {+int(ec1), -int(ec2), 0 }, {+int(ec1), -int(ec2), +int(ec1) }, {+int(ec1), -int(ec2), -int(ec1)}, {+int(ec1), -int(ec2), +int(ec2)}, {+int(ec1), -int(ec2), -int(ec2)},
  {-int(ec1), 0        , 0 }, {-int(ec1), 0        , +int(ec1) }, {-int(ec1), 0        , -int(ec1)}, {-int(ec1), 0        , +int(ec2)}, {-int(ec1), 0        , -int(ec2)},
  {-int(ec1), +int(ec1), 0 }, {-int(ec1), +int(ec1), +int(ec1) }, {-int(ec1), +int(ec1), -int(ec1)}, {-int(ec1), +int(ec1), +int(ec2)}, {-int(ec1), +int(ec1), -int(ec2)},
  {-int(ec1), -int(ec1), 0 }, {-int(ec1), -int(ec1), +int(ec1) }, {-int(ec1), -int(ec1), -int(ec1)}, {-int(ec1), -int(ec1), +int(ec2)}, {-int(ec1), -int(ec1), -int(ec2)},
  {-int(ec1), +int(ec2), 0 }, {-int(ec1), +int(ec2), +int(ec1) }, {-int(ec1), +int(ec2), -int(ec1)}, {-int(ec1), +int(ec2), +int(ec2)}, {-int(ec1), +int(ec2), -int(ec2)},
  {-int(ec1), -int(ec2), 0 }, {-int(ec1), -int(ec2), +int(ec1) }, {-int(ec1), -int(ec2), -int(ec1)}, {-int(ec1), -int(ec2), +int(ec2)}, {-int(ec1), -int(ec2), -int(ec2)},
  {+int(ec2), 0        , 0 }, {+int(ec2), 0        , +int(ec1) }, {+int(ec2), 0        , -int(ec1)}, {+int(ec2), 0        , +int(ec2)}, {+int(ec2), 0        , -int(ec2)},
  {+int(ec2), +int(ec1), 0 }, {+int(ec2), +int(ec1), +int(ec1) }, {+int(ec2), +int(ec1), -int(ec1)}, {+int(ec2), +int(ec1), +int(ec2)}, {+int(ec2), +int(ec1), -int(ec2)},
  {+int(ec2), -int(ec1), 0 }, {+int(ec2), -int(ec1), +int(ec1) }, {+int(ec2), -int(ec1), -int(ec1)}, {+int(ec2), -int(ec1), +int(ec2)}, {+int(ec2), -int(ec1), -int(ec2)},
  {+int(ec2), +int(ec2), 0 }, {+int(ec2), +int(ec2), +int(ec1) }, {+int(ec2), +int(ec2), -int(ec1)}, {+int(ec2), +int(ec2), +int(ec2)}, {+int(ec2), +int(ec2), -int(ec2)},
  {+int(ec2), -int(ec2), 0 }, {+int(ec2), -int(ec2), +int(ec1) }, {+int(ec2), -int(ec2), -int(ec1)}, {+int(ec2), -int(ec2), +int(ec2)}, {+int(ec2), -int(ec2), -int(ec2)},
  {-int(ec2), 0        , 0 }, {-int(ec2), 0        , +int(ec1) }, {-int(ec2), 0        , -int(ec1)}, {-int(ec2), 0        , +int(ec2)}, {-int(ec2), 0        , -int(ec2)},
  {-int(ec2), +int(ec1), 0 }, {-int(ec2), +int(ec1), +int(ec1) }, {-int(ec2), +int(ec1), -int(ec1)}, {-int(ec2), +int(ec1), +int(ec2)}, {-int(ec2), +int(ec1), -int(ec2)},
  {-int(ec2), -int(ec1), 0 }, {-int(ec2), -int(ec1), +int(ec1) }, {-int(ec2), -int(ec1), -int(ec1)}, {-int(ec2), -int(ec1), +int(ec2)}, {-int(ec2), -int(ec1), -int(ec2)},
  {-int(ec2), +int(ec2), 0 }, {-int(ec2), +int(ec2), +int(ec1) }, {-int(ec2), +int(ec2), -int(ec1)}, {-int(ec2), +int(ec2), +int(ec2)}, {-int(ec2), +int(ec2), -int(ec2)},
  {-int(ec2), -int(ec2), 0 }, {-int(ec2), -int(ec2), +int(ec1) }, {-int(ec2), -int(ec2), -int(ec1)}, {-int(ec2), -int(ec2), +int(ec2)}, {-int(ec2), -int(ec2), -int(ec2)},
};
const ftype3 ef_host[Qn] = {
  {0   , 0   , 0 }, {0   , 0   , +ec1 }, {0   , 0   , -ec1}, {0   , 0   , +ec2}, {0   , 0   , -ec2}, // 0--4
  {0   , +ec1, 0 }, {0   , +ec1, +ec1 }, {0   , +ec1, -ec1}, {0   , +ec1, +ec2}, {0   , +ec1, -ec2}, // 5--9
  {0   , -ec1, 0 }, {0   , -ec1, +ec1 }, {0   , -ec1, -ec1}, {0   , -ec1, +ec2}, {0   , -ec1, -ec2}, // 10--14
  {0   , +ec2, 0 }, {0   , +ec2, +ec1 }, {0   , +ec2, -ec1}, {0   , +ec2, +ec2}, {0   , +ec2, -ec2}, // 15--19
  {0   , -ec2, 0 }, {0   , -ec2, +ec1 }, {0   , -ec2, -ec1}, {0   , -ec2, +ec2}, {0   , -ec2, -ec2}, // 20--24
  {+ec1, 0   , 0 }, {+ec1, 0   , +ec1 }, {+ec1, 0   , -ec1}, {+ec1, 0   , +ec2}, {+ec1, 0   , -ec2}, // 25--29
  {+ec1, +ec1, 0 }, {+ec1, +ec1, +ec1 }, {+ec1, +ec1, -ec1}, {+ec1, +ec1, +ec2}, {+ec1, +ec1, -ec2}, // 30--34
  {+ec1, -ec1, 0 }, {+ec1, -ec1, +ec1 }, {+ec1, -ec1, -ec1}, {+ec1, -ec1, +ec2}, {+ec1, -ec1, -ec2}, // 35--39
  {+ec1, +ec2, 0 }, {+ec1, +ec2, +ec1 }, {+ec1, +ec2, -ec1}, {+ec1, +ec2, +ec2}, {+ec1, +ec2, -ec2}, // 40--44
  {+ec1, -ec2, 0 }, {+ec1, -ec2, +ec1 }, {+ec1, -ec2, -ec1}, {+ec1, -ec2, +ec2}, {+ec1, -ec2, -ec2}, // 45--49
  {-ec1, 0   , 0 }, {-ec1, 0   , +ec1 }, {-ec1, 0   , -ec1}, {-ec1, 0   , +ec2}, {-ec1, 0   , -ec2}, // 50--54
  {-ec1, +ec1, 0 }, {-ec1, +ec1, +ec1 }, {-ec1, +ec1, -ec1}, {-ec1, +ec1, +ec2}, {-ec1, +ec1, -ec2}, // 55--59
  {-ec1, -ec1, 0 }, {-ec1, -ec1, +ec1 }, {-ec1, -ec1, -ec1}, {-ec1, -ec1, +ec2}, {-ec1, -ec1, -ec2}, // 60--64
  {-ec1, +ec2, 0 }, {-ec1, +ec2, +ec1 }, {-ec1, +ec2, -ec1}, {-ec1, +ec2, +ec2}, {-ec1, +ec2, -ec2}, // 65--69
  {-ec1, -ec2, 0 }, {-ec1, -ec2, +ec1 }, {-ec1, -ec2, -ec1}, {-ec1, -ec2, +ec2}, {-ec1, -ec2, -ec2}, // 70--74
  {+ec2, 0   , 0 }, {+ec2, 0   , +ec1 }, {+ec2, 0   , -ec1}, {+ec2, 0   , +ec2}, {+ec2, 0   , -ec2}, // 75--79
  {+ec2, +ec1, 0 }, {+ec2, +ec1, +ec1 }, {+ec2, +ec1, -ec1}, {+ec2, +ec1, +ec2}, {+ec2, +ec1, -ec2}, // 80--84
  {+ec2, -ec1, 0 }, {+ec2, -ec1, +ec1 }, {+ec2, -ec1, -ec1}, {+ec2, -ec1, +ec2}, {+ec2, -ec1, -ec2}, // 85--89
  {+ec2, +ec2, 0 }, {+ec2, +ec2, +ec1 }, {+ec2, +ec2, -ec1}, {+ec2, +ec2, +ec2}, {+ec2, +ec2, -ec2}, // 90--94
  {+ec2, -ec2, 0 }, {+ec2, -ec2, +ec1 }, {+ec2, -ec2, -ec1}, {+ec2, -ec2, +ec2}, {+ec2, -ec2, -ec2}, // 95--99
  {-ec2, 0   , 0 }, {-ec2, 0   , +ec1 }, {-ec2, 0   , -ec1}, {-ec2, 0   , +ec2}, {-ec2, 0   , -ec2}, // 100--104
  {-ec2, +ec1, 0 }, {-ec2, +ec1, +ec1 }, {-ec2, +ec1, -ec1}, {-ec2, +ec1, +ec2}, {-ec2, +ec1, -ec2}, // 105--109
  {-ec2, -ec1, 0 }, {-ec2, -ec1, +ec1 }, {-ec2, -ec1, -ec1}, {-ec2, -ec1, +ec2}, {-ec2, -ec1, -ec2}, // 110--114
  {-ec2, +ec2, 0 }, {-ec2, +ec2, +ec1 }, {-ec2, +ec2, -ec1}, {-ec2, +ec2, +ec2}, {-ec2, +ec2, -ec2}, // 115--119
  {-ec2, -ec2, 0 }, {-ec2, -ec2, +ec1 }, {-ec2, -ec2, -ec1}, {-ec2, -ec2, +ec2}, {-ec2, -ec2, -ec2}, // 120--124
};
const ftype W0=(3*TLat*TLat-(ec1d+ec2d)*TLat+ec1d*ec2d)/(ec1d*ec2d);
const ftype W1=TLat*(ec2d-3*TLat)/(2*ec1d*(ec2d-ec1d));
const ftype W2=TLat*(ec1d-3*TLat)/(2*ec2d*(ec1d-ec2d));
const ftype w_host[Qn] = {
  W0*W0*W0, W0*W0*W1, W0*W0*W1, W0*W0*W2, W0*W0*W2,
  W0*W1*W0, W0*W1*W1, W0*W1*W1, W0*W1*W2, W0*W1*W2,
  W0*W1*W0, W0*W1*W1, W0*W1*W1, W0*W1*W2, W0*W1*W2,
  W0*W2*W0, W0*W2*W1, W0*W2*W1, W0*W2*W2, W0*W2*W2,
  W0*W2*W0, W0*W2*W1, W0*W2*W1, W0*W2*W2, W0*W2*W2,
  W1*W0*W0, W1*W0*W1, W1*W0*W1, W1*W0*W2, W1*W0*W2,
  W1*W1*W0, W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W1*W0, W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W2*W0, W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W1*W2*W0, W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W1*W0*W0, W1*W0*W1, W1*W0*W1, W1*W0*W2, W1*W0*W2,
  W1*W1*W0, W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W1*W0, W1*W1*W1, W1*W1*W1, W1*W1*W2, W1*W1*W2,
  W1*W2*W0, W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W1*W2*W0, W1*W2*W1, W1*W2*W1, W1*W2*W2, W1*W2*W2,
  W2*W0*W0, W2*W0*W1, W2*W0*W1, W2*W0*W2, W2*W0*W2,
  W2*W1*W0, W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W1*W0, W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W2*W0, W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
  W2*W2*W0, W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
  W2*W0*W0, W2*W0*W1, W2*W0*W1, W2*W0*W2, W2*W0*W2,
  W2*W1*W0, W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W1*W0, W2*W1*W1, W2*W1*W1, W2*W1*W2, W2*W1*W2,
  W2*W2*W0, W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
  W2*W2*W0, W2*W2*W1, W2*W2*W1, W2*W2*W2, W2*W2*W2,
};
const int reverseX_host[Qn] = {
  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
  25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124
};
const int reverseY_host[Qn] = {
  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
  25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124
};
const int reverseZ_host[Qn] = {
  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
  25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124
};
const int reverseXYZ_host[Qn] = {
  0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
  25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124
};

#elif defined D1Q7
const int Qn=7;

const ftype TLat=0.69795332201968308823840905538932529;///100.;//1/3.;
const ftype Dim=1;

const int3 e_host[Qn] = {
 /*0 -0 */ { 0, 0, 0}, { 1, 0, 0}, { -1, 0, 0}, {2, 0, 0}, {-2, 0, 0}, {3, 0, 0}, {-3, 0, 0}
};
const ftype W0=(36-49*TLat+42*TLat*TLat-15*TLat*TLat*TLat)/36.;
const ftype W1=TLat*(12-13*TLat+5*TLat*TLat)/16.;
const ftype W2=TLat*(-3+10*TLat-5*TLat*TLat)/40.;
const ftype W3=TLat*(4-15*TLat+15*TLat*TLat)/720.;
const ftype w_host[Qn] = {W0,W1,W1,W2,W2,W3,W3};
const int reverseX_host[Qn] = {0, 2, 1, 4,3,6,5};
const int reverseY_host[Qn] = {0, 1, 2,3,4,5,6};
const int reverseZ_host[Qn] = {0, 1, 2,3,4,5,6};
const int reverseXYZ_host[Qn]={0, 2, 1,4,3,6,5};
#endif

//const long int No=(Qn+1)/2;
const long int No=1;

extern __constant__ int3 e[Qn];
extern __constant__ ftype3 ef[Qn];
extern __constant__ ftype w[Qn];
extern __constant__ int reverseX[Qn];  
extern __constant__ int reverseY[Qn];  
extern __constant__ int reverseZ[Qn];  
extern __constant__ int reverseXYZ[Qn];

const ftype cs2 = dx*dx/3.;
const ftype dcs2 = 1./cs2;
const ftype dcs4 = 1./(cs2*cs2);

enum {IndG=0, IndF=1, IndB=2, IndI=3, IndU=4, IndIf=5, IndIg=6, IndS=7};

const int Nxtile = 16;
const int Nytile = 8;

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

#ifdef __CUDA_ARCH__
#define e_c e
#define ef_c ef
#define w_c w
#else
#define e_c e_host
#define ef_c ef_host
#define w_c w_host
#endif
