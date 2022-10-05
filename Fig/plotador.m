C=load('Curves.txt');
nh=314;
figure(1);
for ii=0:287
     plot(C(ii*nh+1:(ii+1)*nh ));
     hold on;
 endfor
 
 
 