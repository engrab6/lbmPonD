from math import *
r1= 1.0
u1= 0.0
us= -2.5#0.5
p1= 1.0
#p1= 3.0
g=5.0/3.0
g=3.0

e1=p1/r1/(g-1)

c1=sqrt(g*p1/r1)

p2 = (((us-u1)**2/c1**2-1)*2*g/(g+1) + 1)*p1
r2 = r1/(1+(p1-p2)/r1/(u1-us)**2)
u2 = ( (r2-r1)*us + r1*u1 )/r2
e2=p2/r2/(g-1)

#print us*(r2-r1), r2*u2-r1*u1
#print us*(r2*u2-r1*u1), r2*u2**2+p2-(r1*u1**2+p1)
#print us*(r2*e2+r2*1/2.*u2**2-r1*e1-r1*1/2.*u1**2), r2*u2*(e2+0.5*u2**2+p2/r2)-r1*u1*(e1+0.5*u1**2+p1/r1)

print r2, u2, p2
print [r1, u1, p1/r1], [r2,u2,p2/r2]

