function y = mult_mov_load_guess (xi, x1, x2, mv, A, freq, x0, t);

i = length(mv);
mv1 = mv(1); A1 = A(1); freq1 = freq(1); x01 = x0(1); 
x = (x1 + x2)/2 + (x2 - x1)*xi/2;
 
b = 0.5;  % Width (standard deviation) of the gauss pulse loading. 

arg = x - x01 - mv1 * t;
spatial_variation  = (cos(arg)+1).*exp(-1*abs(1/b*(arg).^2)); 
 
y1 = -A1*sin(2*pi*freq1*t)   * spatial_variation;   %function of x, t
 

%multiple load
y2 = 0;
if i == 2
mv2 = mv(2); A2 = A(2); freq2 = freq(2); x02 = x0(2);     

arg = x - x02 + mv2 * t;  %CHANGE DIRECTION
 
y2 = -A2*sin(2*pi*freq2*t)   * spatial_variation;   %function of x, t
end

y3 = 0;
if i == 3
 mv3 = mv(3); A3 = A(3); freq3 = freq(3); x03 = x0(3);     

arg = x - x03 - mv3 * t;
 
y3 = -A3*sin(2*pi*freq3*t)   * spatial_variation;   %function of x, t
end


 y = y1+y2+y3; 
 
return
    