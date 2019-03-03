import numpy as np
import itertools
from math import *
#c = list(itertools.product((0,1,-1),(0,1,-1),(0,1,-1)))
c=((0,),(1,),(-1,),(3,),(-3,))
print c
Qn=len(c)
Tl=1.0-sqrt(2./5.)
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

c1=c[1][0]
c2=c[3][0]
  
def gauge_convert(f, ufrom,Tfrom, uto,Tto):
  fnew=np.zeros(Qn)
  e_to = map(lambda cc: cc[0]*sqrt(Tto/Tl), c)
  e_fr = map(lambda cc: cc[0]*sqrt(Tfrom/Tl), c)
  e1_to=e_to[1]
  e2_to=e_to[3]
  e1_fr=e_fr[1]
  e2_fr=e_fr[3]
  for inew in range(len(c)):
    wk = 1./(e1_to**2*e2_to**2) if c[inew][0]==0 else 1./(2*e1_to**2*(e1_to**2-e2_to**2)) if abs(c[inew][0])==c1 else 1./(2*e2_to**2*(e2_to**2-e1_to**2))
    Moment4=33.
    fnew[inew] = Moment4*wk;
    for iold in range(len(c)):
      K4=1
      K3 = -4*e_fr[iold]-e_to[inew]
      if inew==0:
        K2 = 6*e_fr[iold]**2-(e1_to**2+e2_to**2)
        K1 = 2*e_fr[iold]*(e1_to**2+e2_to**2-2*e_fr[iold]**2)
        K0 = (e1_to**2-e_fr[iold]**2)*(e2_to**2-e_fr[iold]**2)
      else:
        K2 = 6*e_fr[iold]**2+3*e_fr[iold]*e_to[inew]-(e1_to**2+e2_to**2-e_to[inew]**2)
        K1 = (2*e_fr[iold]+e_to[inew])*(e1_to**2+e2_to**2-e_to[inew]**2)-e_fr[iold]**2*(4*e_fr[iold]+3*e_to[inew])
        K0 = e_fr[iold]*(e_fr[iold]**2-(e1_to**2+e2_to**2)+e_to[inew]**2)*(e_fr[iold]+e_to[inew])
      #if inew==0:
        #K2 = -(e1_to**2+e2_to**2)
        #K1 = 2*e_fr[iold]*(e1_to**2+e2_to**2)
        #K0 = e1_to**2*e2_to**2-e_fr[iold]**2*e2_to**2-e_fr[iold]**2*e1_to**2
      #else:
        #K2 = 3*e_fr[iold]*e_to[inew]-(e1_to**2+e2_to**2-e_to[inew]**2)
        #K1 = (2*e_fr[iold]+e_to[inew])*(e1_to**2+e2_to**2-e_to[inew]**2)-e_fr[iold]**2*3*e_to[inew]
        #K0 = e_fr[iold]**2*(e_to[inew]**2-(e1_to**2+e2_to**2)) + e_fr[iold]*e_to[inew]*(e_fr[iold]**2+e_to[inew]**2-(e1_to**2+e2_to**2))
      difM4moment = (ufrom[0]+e_fr[iold])**4;
      #M4moment=0
      gx = K4*(uto[0]-ufrom[0])**4 + K3*(uto[0]-ufrom[0])**3 + K2*(uto[0]-ufrom[0])**2 + K1*(uto[0]-ufrom[0]) + K0 - difM4moment
      fnew[inew]+= wk*gx*f[iold]
  return fnew

#f0=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)
f0=(1,2,3,4,5)
f0=(1.273293806252156379,0.362829175487371536,0.362829175487371536,0.000523921386550287,0.000523921386550287)
u0=(0.0,)
T0=1.0

u1=(-0.2,)
T1=1.5

print Moments(u0,T0,f0)
f1 = gauge_convert(f0, u0,T0, u1,T1)
#print f1
print Moments(u1,T1,f1)
print f1

