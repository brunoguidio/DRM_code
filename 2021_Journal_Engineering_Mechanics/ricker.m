function y = ricker(fc,A,t)
%Ricker
     %fc = 10;
     wc = 2*pi*fc;
     tbar = 6*sqrt(6)/wc;
     
     if t <= tbar
         eta = wc*t-3*sqrt(6);
         y = -A*((0.25*eta^2-0.5).*exp(-0.25*eta.^2)-13.*exp(-13.5))...
             ./(0.5+13.*exp(-13.5));
     else
         y = 0;
     end
end