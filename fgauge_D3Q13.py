import numpy as np
from numpy.linalg import inv
import itertools
from math import *
c = list(itertools.product((0,1,-1),(0,1,-1),(0,1,-1)))
_r = sqrt((5+sqrt(5))*0.5)
_s = sqrt((5-sqrt(5))*0.5)
c = (
(0,0,0),
(+_r,+_s,0),(0,+_r,+_s),(+_s,0,+_r),
(+_r,-_s,0),(0,+_r,-_s),(-_s,0,+_r),
(-_r,+_s,0),(0,-_r,+_s),(+_s,0,-_r),
(-_r,-_s,0),(0,-_r,-_s),(-_s,0,-_r),
)
#c=((0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,1),(1,-1),(-1,-1))
print c
Qn=13
#Qn=9
Tl=1.


Momcomb = [(0,0,0),
           (1,0,0),(0,1,0),(0,0,1),
           (2,0,0),(0,2,0),(0,0,2),
           (1,1,0),(1,0,1),(0,1,1),
           #(1,2,0),(1,0,2),(3,0,0),
           (3,0,0),(0,3,0),(0,0,3),
           #(2,2,1),(2,1,2),(1,2,2),
           #(2,1,0),(2,0,1),(1,2,0),(0,2,1),(1,0,2),(0,1,2),
           #(1,1,1),
           #(2,2,0),(2,0,2),(0,2,2),
           #(2,1,1),(1,2,1),(1,1,2),
           #(2,2,2)
           ]
Momcomb_test = [(0,0,0),
           (1,0,0),(0,1,0),(0,0,1),
           #(2,0,0),(0,2,0),(0,0,2),
           #(1,1,0),(1,0,1),(0,1,1),
           #(1,2,0),(1,0,2),(3,0,0),
           #(3,0,0),(0,3,0),(0,0,3),
           #(2,2,1),(2,1,2),(1,2,2),
           #(2,1,0),(2,0,1),(1,2,0),(0,2,1),(1,0,2),(0,1,2),
           #(1,1,1),
           #(2,2,0),(2,0,2),(0,2,2),
           #(2,1,1),(1,2,1),(1,1,2),
           #(2,2,2)
           ]
def Moments(u,T,f):
  retM = np.zeros(Qn)
  imom=0
  if len(Momcomb[0])==2:
    for m,n in Momcomb:
      for i in range(len(c)):
        retMi = (sqrt(T/Tl)*c[i][0]+u[0])**m*(sqrt(T/Tl)*c[i][1]+u[1])**n
        retM[imom]+= retMi*f[i]
      imom+=1
  elif len(Momcomb_test[0])==3:
    for m,n,p in Momcomb_test:
      for i in range(len(c)):
        retMi = (sqrt(T/Tl)*c[i][0]+u[0])**m*(sqrt(T/Tl)*c[i][1]+u[1])**n*(sqrt(T/Tl)*c[i][2]+u[2])**p
        retM[imom]+= retMi*f[i]
      imom+=1
    for iTflux in (0,1,2):
      for i in range(len(c)):
        Tmpr = (sqrt(T/Tl)*c[i][0]+u[0])**2+(sqrt(T/Tl)*c[i][1]+u[1])**2+(sqrt(T/Tl)*c[i][2]+u[2])**2
        retMi = (sqrt(T/Tl)*c[i][iTflux]+u[iTflux])*Tmpr;
        retM[imom]+= retMi*f[i]
      imom+=1
  return retM

def gauge_convert_inv(f, ufrom,Tfrom, uto,Tto):
  fnew=np.zeros(Qn)
  Mto   = np.zeros(shape=(Qn,Qn))
  Mfrom = np.zeros(shape=(Qn,Qn))
  imom=0
  if len(c[0])==2:
    for m,n in Momcomb:
      for i in range(Qn):
        Mfrom[imom][i] = (sqrt(Tfrom/Tl)*c[i][0]+ufrom[0])**m*(sqrt(Tfrom/Tl)*c[i][1]+ufrom[1])**n
        Mto[imom][i]   = (sqrt(Tto/Tl)*c[i][0]+uto[0])**m*(sqrt(Tto/Tl)*c[i][1]+uto[1])**n
      imom+=1
  elif len(c[0])==3:
    for m,n,p in Momcomb:
      for i in range(Qn):
        Mfrom[imom][i] = (sqrt(Tfrom/Tl)*c[i][0]+ufrom[0])**m*(sqrt(Tfrom/Tl)*c[i][1]+ufrom[1])**n*(sqrt(Tfrom/Tl)*c[i][2]+ufrom[2])**p
        Mto[imom][i]   = (sqrt(Tto/Tl)*c[i][0]+uto[0])**m*(sqrt(Tto/Tl)*c[i][1]+uto[1])**n*(sqrt(Tto/Tl)*c[i][2]+uto[2])**p
      imom+=1
    #for iTflux in (0,1,2):
      #for i in range(len(c)):
        #fromTmpr = (sqrt(Tfrom/Tl)*c[i][0]+ufrom[0])**2+(sqrt(Tfrom/Tl)*c[i][1]+ufrom[1])**2+(sqrt(Tfrom/Tl)*c[i][2]+ufrom[2])**2
        #toTmpr   = (sqrt(Tto  /Tl)*c[i][0]+uto  [0])**2+(sqrt(Tto  /Tl)*c[i][1]+uto  [1])**2+(sqrt(Tto  /Tl)*c[i][2]+uto  [2])**2
        #Mfrom[imom][i]  = (sqrt(Tfrom/Tl)*c[i][iTflux]+ufrom[iTflux])*fromTmpr;
        #Mto  [imom][i]  = (sqrt(Tto  /Tl)*c[i][iTflux]+uto  [iTflux])*toTmpr;
      #imom+=1      
  #print Mto
  Minv = inv(Mto);
  #print np.dot(Mto, Minv) 
  print np.linalg.det(Mto)
  #print Minv
  fnew = np.dot( Minv, np.dot(Mfrom, f) )
  return fnew

def gauge_convert(f, ufrom,Tfrom, uto,Tto):
  fnew=np.zeros(Qn)
  for inew in range(len(c)):
    wk = -1./Tto if c[inew][0]==0 else 1./(2*Tto)
    wl = -1./Tto if c[inew][1]==0 else 1./(2*Tto)
    wm = -1./Tto if c[inew][2]==0 else 1./(2*Tto)
    for iold in range(len(c)):
      Ax_i = (uto[0]-ufrom[0])/sqrt(3.)-c[iold][0]*sqrt(Tfrom)
      Ay_j = (uto[1]-ufrom[1])/sqrt(3.)-c[iold][1]*sqrt(Tfrom)
      Az_p = (uto[2]-ufrom[2])/sqrt(3.)-c[iold][2]*sqrt(Tfrom)
      Bx_ik = Tto if c[inew][0]==0 else c[inew][0]*sqrt(Tto)*Ax_i
      By_jl = Tto if c[inew][1]==0 else c[inew][1]*sqrt(Tto)*Ay_j
      Bz_pm = Tto if c[inew][2]==0 else c[inew][2]*sqrt(Tto)*Az_p
      gx=Ax_i**2-Bx_ik
      gy=Ay_j**2-By_jl
      gz=Az_p**2-Bz_pm
      fnew[inew]+= wk*wl*wm*gx*gy*gz*f[iold]
  return fnew

f0=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)
f0=(1,2,3,4,5,6,7,8,9,10,11,12,13)
u0=(1.,2.,3.)
T0=3.

u1=(0.,0.,0.)
T1=1./3.

print Moments(u0,T0,f0)
f1 = gauge_convert_inv(f0, u0,T0, u1,T1)
print f1
print Moments(u1,T1,f1)

