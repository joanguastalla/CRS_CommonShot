import numpy as np
from crs_commonshot import *
import matplotlib.pyplot as plt
import time
#import multiprocessing as mp
with open('/home/joan/Pluto_Dataset/pluto.bin','rb') as f:
	data=f.read()
	U=np.frombuffer(data,dtype='f')


# Data sizes(f: Ricker Wavelet peak frequency in Hz)
ns=1126
nh=338
nshots=694
f=15

# Function crs_cs required parameters
U=np.reshape(U,[ns,nh*nshots],order='F')
ds=150
dh=75/2
h=(12*dh + np.arange(0,338*dh,dh))/2
s=np.arange(0,694*ds,ds)
w=np.sqrt(6)/(np.pi*f)
dt=0.008
t=np.arange(0,1126*dt,dt)
alpha=20
dalpha=1
v0=1.5
vmin=-1.5
vmax=1.8
hyperbola_jump=2

# Scale transform feet to km unit
feet2m=.3048
m2km=0.001
s,h=s*feet2m*m2km,h*feet2m*m2km
ds,dh=ds*feet2m*m2km,dh*feet2m*m2km


Semblance=np.zeros([ns*nshots,])
ntimes=1
start=time.time()
for jj in [400]:
    Semblance[jj]=crs_cs(s[200],t[jj],U,h,s,ds,dh,w,dt,alpha,dalpha,v0,vmin,vmax,hyperbola_jump)
end=time.time()
print("Elapsed time is: {} seconds".format(end-start))

