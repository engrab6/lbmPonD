import numpy as np
import itertools
from math import *
#c = list(itertools.product((0,1,-1),(0,1,-1),(0,1,-1)))
#c=((0,0),(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,1),(1,-1),(-1,-1))
#c=((0,),(1,),(-1,),(2,),(-2,),(3,),(-3,))
c=((0,),(2,),(-2,))
c=((0,),(1,),(-1,),(3,),(-3,))
#c=((0,),(1,),(-1,))
print c
Qn=len(c)
Tl=4./3
#Tl=1./3.*(2-(7./(25+15*sqrt(30)))**(1./3)+(5./7.+3./7.*sqrt(30))**(1./3)/5**(2./3.))
print Tl
def Moments(u,T,f):
  retM = np.zeros(Qn*Qn)
  imom=0
  for m in (0,1,2,3,4,5,6):
    for i in range(len(c)):
      retMi = (sqrt(T/Tl)*c[i][0]+u[0])**m
      retM[imom]+= retMi*f[i]
    imom+=1
  return retM

def gauge_convert(f, ufrom,Tfrom, uto,Tto):
  fnew=np.zeros(Qn)
  for inew in range(len(c)):
    wk = -1./Tto if c[inew][0]==0 else 1./(2*Tto)
    #wl = -1./Tto if c[inew][1]==0 else 1./(2*Tto)
    #wm = -1./Tto if c[inew][2]==0 else 1./(2*Tto)
    for iold in range(len(c)):
      Ax_i = (uto[0]-ufrom[0])*sqrt(Tl)-c[iold][0]*sqrt(Tfrom)
      #Ay_j = (uto[1]-ufrom[1])*sqrt(Tl)-c[iold][1]*sqrt(Tfrom)
      #Az_p = (uto[2]-ufrom[2])*sqrt(Tl)-c[iold][2]*sqrt(Tfrom)
      Bx_ik = Tto if c[inew][0]==0 else 1./c[inew][0]*sqrt(Tto)*Ax_i
      #By_jl = Tto if c[inew][1]==0 else c[inew][1]*sqrt(Tto)*Ay_j
      #Bz_pm = Tto if c[inew][2]==0 else c[inew][2]*sqrt(Tto)*Az_p
      gx=((uto[0]-ufrom[0])*sqrt(Tl)-c[iold][0]*sqrt(Tfrom))**2-Bx_ik
      if(c[iold][0]!=0): gx = ((uto[0]-ufrom[0])*sqrt(Tl)-c[iold][0]*sqrt(Tfrom))**2/4.-Bx_ik
      #gy=Ay_j**2-By_jl
      #gz=Az_p**2-Bz_pm
      fnew[inew]+= wk*gx*f[iold]
      #fnew[inew]+= wk*wl*wm*gx*gy*gz*f[iold]
  return fnew

#f0=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)
f0=(1,2,3,4,5,6,7)
u0=(0.,)
T0=1.

u1=(1.,)
T1=3.

print Moments(u0,T0,f0)
f1 = gauge_convert(f0, u0,T0, u1,T1)
print f1
print Moments(u1,T1,f1)

