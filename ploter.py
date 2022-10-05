import numpy as np
from seisplot import seisplot
import matplotlib.pyplot as plt
import time
#import multiprocessing as mp
with open('/home/joan/C_scripts/ProjetoCRS_C/Fig/Pluto.txt','r') as f:
	U=[[float(num) for num in line.split(' ')] for line in f]
print(U[0,0:10]);
#ns=1126
#nh=314
#nshots=694
#U=np.reshape(U,[ns,nh*nshots],order='F')
