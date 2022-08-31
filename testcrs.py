import numpy as np
from crs_commonshot import *
with open('/home/joan/Pluto_Dataset/pluto.bin','rb') as f:
	data=f.read()
	U=np.frombuffer(data,dtype='f')
ns=1126
ds=150
nh=338
dh=75
nshots=694
dt=0.008
v0=1500
s=np.arange(0,694*ds,ds)
h=12*dh + np.arange(0,338*dh,dh)
vmin=1500
vmax=4000
f=15
alpha=30
dalpha=0.1
w=np.sqrt(6)/(np.pi*f)
U=np.reshape(U,[ns,nh*nshots],order='F')
Semblance=np.zeros([ns,nshots])
time=np.arange(0,ns*dt,dt)
meter2feet=3.28084
# Velocities in feet/s
vmin,vmax,v0=vmin*meter2feet,vmax*meter2feet,v0*meter2feet
for ii in range(len(s)):
	for jj in range(ns): 
		Semblance[ii,jj]=CRSCS(U,h,s,w,s[ii],time[jj],dt,alpha,dalpha,v0,vmin,vmax)

