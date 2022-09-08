import numpy as np
from math import isinf
from cubicinterpol import *
from common_receiver import common_receiver


def crs_cs(m0,t0,U,h,s,ds,dh,w,dt,alphamin,alphamax,dalpha,v0,tdown,tup,hyperbola_jump):
    """
    m0        -  Central zero-offset ray midpoint
    t0        -  Time of propagation of the central zero-offset ray
    U         -  All Common shot sections
    h         -  Vector of offsets 
    s         -  Vector of shots
    w         -  Wavelet duration(s)
    dt        -  Time Sampling(s)
    alphamin  -  Minimum emergence angle in degrees of ZO trace (m0,t0)  
    alphamax  -  Maximum emergence angle in degrees of ZO trace (m0,t0)
    dalpha    -  Angle sampling for hyperbola fitting
    v0        -  Velocity at ZO location
    tdown     -  Minimum time, associated with hyperbola fitting  at h[end], given by (t0-tdown)
    tup       -  Maximum time, associated with hyperbola fitting  at h[end], given by (t0+tup)  
    hyperbola_jump  - Number of samples in diference of time in the hyperbola fitting
                      between h[end] and h[0] 
                   

    Computes DCRS parameters, associated to a NIP diffraction point, from a given zero-offset
    trace and traveltime (m0,t0). Configuration of the traces is supossed to be Common Shot,
    linearly distributed with traces sorted as:

    U[:,index]= trace parametrized by the functions s(index) and h(index), calculated by

    s(index)=(index//len(h))*ds
    h(index)=h[0] + mod(index,len(h))*dh

    where,

    ds=s[1]-s[0]
    dh=h[1]-h[0]

    """
    ndown=min(int(t0/dt),int(tdown/dt))
    nup=int(tup/dt)
    nsamples=nup+ndown
    emerge=np.radians([alphamin,alphamax,dalpha])
    shotind=int((m0/ds)*len(h))
    stack_traces=np.arange(shotind,shotind+len(h))
    
# Reciprocity principle gives negative offsets for CS with s=m0
    
    h,geos=common_receiver(2*h,s,m0)
    h_s=np.power(h,2)
    stack_traces=np.insert(stack_traces,0,geos)
   
# Windowing around diffraction traveltime 



    nsw=w/dt
    nsw_lim=nsw//2
    janela=np.arange(-nsw_lim+1,nsw_lim+1)
    factor=2/v0
    vel=np.zeros([2,])
    vel_max=np.zeros([2,])
    alfa_max=np.zeros([2,])
    semblance_max=0
    trsum=np.empty([len(janela),],dtype='float')
    trsum2=0

# Best fitting hyperbola parameters (alfa,vel), within CS section s=m0 

### Compute alfa,vel,semblance by means of length 2 vectors

## Inside while loop
#    vel[0] - always max semblance associated velocity before vel[1]
#    vel[1] - actual iterating velocity

## Inside for loop
#    Index 1 - Gives alfa,vel[0],Semblance 
#    Index 0 - Compare semblance_alfa[0] with Semblance
#       if Semblance > semblance_alfa[0]       then: alfa_max[0]=alfa
#                                              vel_max[0]=v[0]             
#                                              semblance_max[0]=Semblance
# mid_aperture=vel[1]/2*np.sqrt(t0*w)
# h_aperture=h[abs[h]<mid_aperture]
  
    for alfa in np.arange(emerge[0],emerge[1]+emerge[2],emerge[2]): 
        scale=factor*np.sin(alfa)
        Semblance=0
        tdifrac=t0*(t0 + 2*scale*h) + scale**2*h_s
        for ii in range(0,nsamples,hyperbola_jump):
            vel[1]=((t0 + (ii-ndown)*dt)**2 - (t0+scale*h[-1])**2)/(h[-1]**2)
            tdifrac=np.sqrt(tdifrac + vel[1]*h_s)
            hyper_stack=np.array([cubic_interpol(U[:,stack_traces],dt,tdifrac + dt*window) for window in janela])
            trsum=np.sum(np.power(np.sum(hyper_stack,axis=1),2))
            trsum2=np.sum(np.power(hyper,2))
            if isinf(np.sum(np.power(trsum,2))):
                    print("happens for vel:{},alfa:{} and t0:{}".format(vel[1],alfa,t0))
                    return -1  
            semb=trsum/(len(h)*trsum2)
            indmax=np.argmax([Semblance,semb])
            Semblance=np.max([Semblance,semb])
            vel[0]=vel[indmax]
        vel_max[1]=vel[0]
        alfa_max[1]=alfa
        semblance_max=max(semblance_max,Semblance)
        indmax_alfa=np.argmax([semblance_max,Semblance])
        vel_max[0]=vel_max[indmax_alfa]
        alfa_max[0]=alfa_max[indmax_alfa]
        

    return semblance_max
