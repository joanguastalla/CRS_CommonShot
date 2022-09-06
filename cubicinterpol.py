import numpy as np
def cubic_interpol(s,dt,time):
	"""
	s -  trace to be interpoled
	dt - trace sampling
	t  - time interpolation vector
	
	Funtion that does cubic interpolation, using as conditions:
	
	s[start],s[end],s'[start],s'[end] 	
	
	where,
	
	start=floor(t/dt), end=ceil(t/dt)
	
    and s'[start],s'[end] are centered finite diffrence derivatives, i.e., 
	
	s'[start]=(s[end]-s[start-1])/(2*dt), s'[end]=(s[end+1]-s[start])/(2*dt) 

	"""	
	interpol=np.empty([0,])
	alpha=np.empty([2,1],dtype=float)
	A=np.zeros([4,4],dtype=float)
	b=np.empty([4,1],dtype=float)
	powers=np.arange(3,-1,-1)
	for t in time:
		tt=t/dt
		ind=max(0,int(np.floor(tt)))
		ind=min(ind,len(s)-2)
		if abs(np.ceil(t/dt))==int(t/dt) or t<0:
			interpol=np.append(interpol,s[ind])
			continue
		else:
			alpha[0]=s[ind+1]-s[ind-1]
			alpha[1]=s[ind+2]-s[ind]
			alpha=alpha/(2*dt)
			A[0,:]=np.power(ind*dt,powers)[np.newaxis]
			A[1,:]=np.power((ind+1)*dt,powers)[np.newaxis]
			A[2,:-1]=powers[0:3]*np.power(ind*dt,powers[1:])[np.newaxis]
			A[-1,:-1]=powers[0:3]*np.power((ind+1)*dt,powers[1:])[np.newaxis]
			b[0]=s[ind]
			b[1]=s[ind+1]
			b[2]=alpha[0]
			b[-1]=alpha[1]
			c=np.linalg.solve(A,b)
			c=np.reshape(c,[len(c),])
			interpol=np.append(interpol,np.sum(c*np.power(t,powers)))
	return interpol

