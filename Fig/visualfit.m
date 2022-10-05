dt=0.008;
feet2km=0.0003048;
dh=75;
dh=dh*feet2km;
nh=100;
h0=36*75*feet2km;
h=h0:dh:(h0+(nh-1)*dh);
hminus=-flip(h);
hfull=[hminus h];
v0=1.5;
data=U(:,(301*314+1):(301*314 + nh));
CRSZO=zeros(100,1);
ep=40;
#l=nh/4 -1; 
l=nh-1;
sided=1;
p=4;
c_stack=hannfilter(l,ep,p,sided);
CRS=load("CRS.txt");
t0=130*dt;
alfa=CRS(1,2);
vel=CRS(1,3);
vel=vel/t0;
parabolic=t0 + sin(alfa)/v0*hfull + vel*hfull.^2;
figure(4);
plot(parabolic);
set(gca,'Ydir','reverse');
for jj=1:1:870
 t0=(129+jj)*dt;
  alfa=CRS(jj,2);
  vel=CRS(jj,3);  
  hypertime=(t0 + sin(alfa)/v0*h).^2 + vel*h.^2;
  hypertime=sqrt(hypertime);
  for kk=1:1:nh
    #aux_stack(kk)=cubicinterpol(data(:,kk),dt,hypertime(kk)); 
      ind=floor(hypertime(kk)/dt);
      aux_stack(kk)=(data(ind,kk)+data(ind+1,kk))/2;
  endfor
  CRSZO(jj)=2*sum(aux_stack.*c_stack)/nh;
 endfor
figure(2)
 plot(CRSZO);
hold on;
plot(ZO(130:999,301));
 #imagesc(data(131:230,:));
figure(1);
imagesc(data(1:500,:));
hold on;
for ii=1:1:100
 if CRS(ii,1)>0.1
  t0=(129+ii)*dt;
  alfa=CRS(ii,2);
  vel=CRS(ii,3);  
  hypertime=(t0 + sin(alfa)/v0*h).^2 + vel*h.^2;
  hypertime=sqrt(hypertime);
  hypertime=round(hypertime/dt);
  plot(hypertime(1:nh));
 endif
endfor