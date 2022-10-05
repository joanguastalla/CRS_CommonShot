## Copyright (C) 2022 Joan Guastalla Arruda

## Author: Joan Guastalla Arruda <joan@lince>
## Created: 2022-09-30

function amp= cubicinterpol (data,dt,t)
ind=floor(t/dt)+1;
if ind==(t/dt +1)
   amp=data(ind);
else
  ind=min(ind,length(data)-1);  
  if ind>1
    alfa0=(data(ind+1)-data(ind-1))/(2*dt);
  else
    alfa0=(data(ind+1)-data(ind))/dt;
  endif
  
  if ind<(length(data)-1) 
    alfa1=(data(ind+2)-data(ind))/(2*dt); 
  else
    alfa1=(data(ind+1)-data(ind))/dt;
  endif
  t0=ind*dt;
  t1=t0+dt;
  A=[[t0^3,t0^2,t0,1];[t1^3,t1^2,t1,1];
       [3*t0^2,2*t0,1,0];[3*t1^2,2*t1,1,0]];
  b=[data(ind),data(ind+1),alfa0,alfa1];
  b=b';
  x=linsolve(A,b);
  amp=x(1)*t^3 + x(2)*t^2 + x(3)*t + x(4);
 endif
endfunction
