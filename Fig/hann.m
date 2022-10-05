function hanning = hann(N)
  for ii=0:1:N
    hanning(ii+1)=0.5*(1-cos(2*pi*ii/N));
   endfor
   hanning=hanning';
 endfunction