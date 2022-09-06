import numpy as np
def common_receiver(h,s,g):
    """
    h - offsets
    s - shots 
    g - geophone

    Given a common shot section, this code computes 
    all possible pair (s,h) such that g=s+2h, i.e., 
    a common receiver section

    Returned parameters are offsets and traces in-
    dexes to get the CR section from the CS data.
    """
    ds=s[1]-s[0]
    dh=h[1]-h[0]
    ratio=ds//dh
    shotind=(g-s[0])//ds
    gind=max(0,shotind - int((len(h)-1)/ratio)-h[0]//ds)
    geos=np.arange(gind,shotind-h[0]//ds)
    geos=geos*len(h) + ratio*np.arange(len(geos)-1,-1,-1)
    h=np.insert(h,0,-np.arange(len(geos)-1,-1,-1)*ratio*dh - h[0])

    
    return h,geos


