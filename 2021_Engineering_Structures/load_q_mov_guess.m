function y = load_q_mov_guess(xi, x1, x2, mv, A, freq, x0, t)

 
x = (x1 + x2)/2 + (x2 - x1)*xi/2;
 
%x0 = 50; % initial position of the central point of loading. 
%mv = 20 ;% moving speed of the loading [m/s]
b = 0.5;  % Width (standard deviation) of the gauss pulse loading.
%A = -100; %Amplitude 
%freq = 50; % [Hz]

arg = x - x0 - mv * t; % relative location of Gauss quadrature point, xi, to the center of the i-th loading
spatial_variation  = (cos(arg)+1).*exp(-1*abs(1/b*(arg).^2)); 
 
 %y = -A*sin(2*pi*freq*t)   * spatial_variation;   %function of x, t
  y = spatial_variation; 
return
    