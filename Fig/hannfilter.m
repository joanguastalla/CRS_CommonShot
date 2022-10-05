function [c_stack]=hannfilter(l,ep,p,sided)
     c_stack=hann(2*ep+1).^p;
     if sided==0
        c_stack=[c_stack(1:ep)' ones(1,2*(l-ep)+1) c_stack(ep+2:end)'];
      else
        c_stack= [ones(1,(l-ep)) c_stack(ep+2:end)'];
      endif
end   