import matplotlib.pyplot as plt
import numpy as np
def seisplot(U,t0=0,dt=0.004,nh=1,ftr=0,dtr=1,ntr=0,f2=0,d2=1,nlabel=8,clip=0.1,interpolation='spline36'):
    """ 
    U      -   Matrix of traces
    t0     -   Initial time of ploting
    dt     -   Time sampling
    nh     -   Number of offsets
    ftr    -   First trace to show
    dtr    -   Sampling of trace ploting
    f2     -   First coordinate 
    d2     -   Coordinate sampling
    nlabel -   Number of label marks(+-1) 
    clip   -   Between (0,1] is the percentile for clipping max amplitude on U
    """
    ind0=int(t0/dt)
    size=np.shape(U)
    if ntr==0:
        ntr=(size[1]-ftr)/dtr + 1
    
    if dtr <= 0 or type(dtr)!=int:
        raise('dtr must be integer greater than 0')
    
    try:
        U=U[ind0:,:]
    except:
        Exception('Initial time t0 must be between 0 and ',(size[0]-1)*dt)
    
    
    try:
        aux=np.arange(int(ftr),int(ntr)*int(dtr),int(dtr))
        U=U[:,aux]
    except:
        Exception('ftr,ntr or dtr not set properly')
    size=np.shape(U)    
    max_value=np.max(np.max(np.abs(U)))*clip
    ax=plt.gca()
    plt.imshow(U,cmap='gist_yarg',vmin=-max_value,vmax=max_value,aspect='auto')
    xticks=np.arange(0,size[1],int(size[1]/nlabel))
    xtick_labels=100*((xticks//nh)*d2)
    xtick_labels=xtick_labels.astype(int)/100
    ax.set_xticks(xticks)
    ax.set_xticklabels(xtick_labels)
    yticks=np.arange(0,size[0],100)
    ytick_labels=10*(yticks*dt)
    ytick_labels=ytick_labels.astype(int)/10
    ax.set_yticks(yticks)
    ax.set_yticklabels(ytick_labels)
    plt.show()
