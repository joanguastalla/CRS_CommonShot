import numpy as np
from cubicinterpol import *
from common_receiver import common_receiver

def CRSCS(U,h,s,w,m0,t0,dt,alpha,dalpha,v0,vmin,vmax):
    """
    U     -  All Common shot sections
    h     -  Vector of offsets 
    s     -  Vector of shots
    w     -  Wavelet duration(s)
    t0    -  Time(s)
    dt    -  Time Sampling(s)
    alpha -  Maximum emergence angle of ZO trace (m0,t0)
    v0    -  Velocity at ZO location
    vmin  -  Minimum curvature of hyperbolic fitting
    vmax  -  Maximum curvature of hyperbolic fitting
    
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

    size=np.shape(U)
    ns=size[0]
    ds=s[1]-s[0]
    dh=h[1]-h[0]
    alpha=alpha*np.pi/180
    dalpha=dalpha*np.pi/180
    shotind=int((m0/ds)*len(h))
    stack_traces=np.arange(shotind,shotind+len(h))
    
# Reciprocity principle gives negative offsets for CS with s=m0
    
    h,geos=common_receiver(h,s,m0)
    
    stack_traces=np.insert(stack_traces,0,geos)

   
# Windowing around diffraction traveltime 



    nsw=w/dt
    nsw_lim=nsw//2
    if np.mod(nsw,2)==0:
        janela=np.arange(-nsw_lim+1,nsw_lim+1)
    else:
        janela=np.arange(-nsw_lim,nsw_lim+1)
    factor=2/v0
    alfas=np.arange(-alpha,alpha+dalpha,dalpha)
    vel=np.array([0,vmin])
    vel_max=np.zeros([2,])
    alfa_max=np.zeros([2,])
    semblance_alfa=np.zeros([2,])
    trsum=np.empty([len(janela),],dtype='float')
    trsum2=0

# Best fitting hyperbola parameters (alfa,vel), within CS section s=m0 

### Compute alfa,vel,semblance by means of length 2 vectors

## Inside while loop
#    vel[0] - always max semblance associated velocity before v[1]
#    vel[1] - actual iterating velocity

## Inside for loop
#    Index 1 - Gives alfa,vel[0],Semblance 
#    Index 0 - Compare semblance_alfa[0] with Semblance
#       if semblance_alfa[0] > Semblance then: alfa_max[0]=alfa
#                                              vel_max[0]=v[0]             
#                                              semblance_max[0]=Semblance
# mid_aperture=vel[1]/2*np.sqrt(t0*w)
# h_aperture=h[abs[h]<mid_aperture]
    for alfa in alfas: 
        scale=factor*np.sin(alfa)
        Semblance=0
        while abs(vel[1]) < vmax:
            for hh in range(len(h)):
                 tdifrac=t0/2 + np.sqrt((t0/2 + scale*h[hh])**2 + \
                 (h[hh]/vel[1])**2)
                 aux=cubic_interpol(U[:,stack_traces[hh]],dt,tdifrac + \
                 dt*janela)
                 trsum+=aux
                 trsum2+=np.sum(np.power(aux,2))
            semb=np.sum(np.power(trsum,2))/(len(h)*trsum2)
            indmax=np.argmax([Semblance,semb])
            Semblance=np.max([Semblance,semb])
            vel[0]=vel[indmax]
            tmax= t0/2 + np.sqrt((t0/2 + scale*h[-1])**2  + \
            (2*h[-1]/vel[1])**2)
            vel[1]=np.sqrt(1/((1/vel[1])**2 - 2*dt*tmax/4*(h[-1]**2) + \
            dt**2/(4*h[-1]**2)))
        vel_max[1]=vel[0]
        alfa_max[1]=alfa
        semblance_alfa[1]=Semblance
        indmax_alfa=np.argmax(semblance_alfa)
        semblance_alfa[0]=semblance_alfa[indmax_alfa]
        vel_max[0]=vel_max[indmax_alfa]
        alfa_max[0]=alfa_max[indmax_alfa]
        

    return semblance_alfa[0],alfa_max[0],vel_max[0]
