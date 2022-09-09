import numpy as np
from crs_commonshot import *
from itertools import repeat
import concurrent.futures as parallel
import matplotlib.pyplot as plt
import time
from seisplot import seisplot
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
dt=0.008
t=np.arange(0,1126*dt,dt)
alphamin=-10
alphamax=20
dalpha=2
v0=1.5
tdown=0
tup=5
hyperbola_jump=40
chunk_best=25
# Scale transform feet to meter unit
feet2km=.0003048
s,h=s*feet2km,h*feet2km/2
ds,dh=ds*feet2km,dh*feet2km


Semblance=np.zeros([ns*nshots,])
def auxcrs(arg):
    crs_cs(s[int(arg//ns)],t[int(np.mod(arg,ns))],U,h,s,ds,dh,w,dt,alphamin,alphamax,dalpha,v0,tdown,tup,hyperbola_jump)

seisplot(U[:,np.arange(100*len(h),101*len(h))],t0=0,dt=0.008)

def main(executor,iterator):
    return executor.map(auxcrs,iterator,chunksize=25)
if __name__=="__main__":
    start_time=time.time()
    ntimes=range(400)
    count=0
    with parallel.ProcessPoolExecutor(max_workers=16) as executor:
        print("Number of workers: ",executor._max_workers)
        Semblance=main(executor,ntimes)
    for semb in Semblance:
        count+=1
    print("Size of semblance: {}".format(count))
    print("Elapsed time: {} seconds".format(time.time()-start_time))
    # Plotting zero-offset traces
        
