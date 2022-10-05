dt=0.008;
nt=1126;
trace_original=load("Interpolated.txt");
for ii=1:1:4*nt-3
  trace_interpoled(ii)=cubicinterpol(trace_original,dt,(ii-1)*dt/4);
endfor
figure(1);
plot(trace_original);
hold on;
plot(trace_interpoled);