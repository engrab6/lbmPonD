import numpy as np
import itertools
from math import *
c = list(itertools.product((0,1,-1),(0,1,-1)))
#c = list(itertools.product((0,1,-1),(0,1,-1),(0,1,-1)))
#c=((0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,1),(1,-1),(-1,-1))
print c
Qn=9
Tl=1./3
def Moments(u,T,f):
  retM = np.zeros(Qn)
  imom=0
  for m,n in (0,0),(1,0),(0,1),(1,1),(2,0),(0,2),(2,1),(1,2),(2,2):
    for i in range(len(c)):
      retMi = (sqrt(T/Tl)*c[i][0]+u[0])**m*(sqrt(T/Tl)*c[i][1]+u[1])**n
      retM[imom]+= retMi*f[i]
    imom+=1
  return retM

def gauge_convert(f, ufrom,Tfrom, uto,Tto):
  fnew=np.zeros(Qn)
  for inew in range(len(c)):
    wk = 1./Tto if c[inew][0]==0 else -1./(2*Tto)
    wl = 1./Tto if c[inew][1]==0 else -1./(2*Tto)
    for iold in range(len(c)):
      Ax_i = (uto[0]-ufrom[0])/sqrt(3.)-c[iold][0]*sqrt(Tfrom)
      Ay_j = (uto[1]-ufrom[1])/sqrt(3.)-c[iold][1]*sqrt(Tfrom)
      Bx_ik = Tto if c[inew][0]==0 else c[inew][0]*sqrt(Tto)*Ax_i
      By_jl = Tto if c[inew][1]==0 else c[inew][1]*sqrt(Tto)*Ay_j
      gx=Ax_i**2-Bx_ik
      gy=Ay_j**2-By_jl
      fnew[inew]+= wk*wl*gx*gy*f[iold]
  return fnew

f0=(1,2,3,4,5,6,7,8,9)
u0=(1.,2.)
T0=2.

u1=(3.,4.)
T1=1.

print Moments(u0,T0,f0)
f1 = gauge_convert(f0, u0,T0, u1,T1)
print f1
print Moments(u1,T1,f1)

