ind=1;
dt=0.008;
t=0.009;
for ii=1:2
    for jj=1:4
      A(ii,jj)=((ind+ii-1)*dt)^(4-jj)
      endfor
endfor

for ii=3:4
    for jj=1:4
      A(ii,jj)=(4-jj)*((ind+ii-3)*dt)^(3-jj)
    endfor
endfor