import numpy as np
from seisplot import seisplot
from crs_commonshot import *
import matplotlib.pyplot as plt
import time
#import multiprocessing as mp
with open('/home/joan/Pluto_Dataset/pluto.bin','rb') as f:
	data=f.read()
	U=np.frombuffer(data,dtype='f')


# Data sizes(f: Ricker Wavelet peak frequency in Hz)
ns=1126
nh=314
nshots=694
f=15

# Function crs_cs required parameters
U=np.reshape(U,[ns,nh*nshots],order='F')
ds=150
dh=75
h=36*dh + np.arange(0,314*dh,dh)
s=np.arange(0,694*ds,ds)
w=np.sqrt(6)/(np.pi*f)
print("Wavelet Time: ",w)
dt=0.008
t=np.arange(0,1126*dt,dt)
alphamin=-10
alphamax=20
dalpha=1
v0=1.5
tdown=0.5
tup=1.0
hyperbola_jump=16

#Parameter for plotting
m0=s[100]
indshot=int(m0/ds)
stack_traces=np.arange(indshot*len(h),indshot*len(h) + len(h))
h,geos=common_receiver(h,s,ds,dh,m0)
stack_traces=np.insert(stack_traces,0,geos)
seisplot(U[:,stack_traces],dt=0.008)
# Scale transform feet to km unit
feet2km=.0003048
s,h=s*feet2km,h*feet2km
ds,dh=ds*feet2km,dh*feet2km


Semblance=np.zeros([ns*nshots,])
ntimes=1
start=time.time()
for jj in [400]:
    Semblance[jj]=crs_cs(s[200],t[jj],U,h,s,ds,dh,w,dt,alphamin,alphamax,dalpha,v0,tdown,tup,hyperbola_jump)
end=time.time()
print("Elapsed time is: {} seconds".format(end-start))

