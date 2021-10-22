function [grad] = compute_grad(F_g,l_bottom) %BG 06/11 - version b
global h; global dt; global R_factor;global R_factor_intensity;

%% Compute the first part
[rows cols] = size(F_g);
for i = 1:rows
    for j = 1:cols
    %finite difference
        if (j == 1)
        %forward difference
        d2Fdx2(i,j) = (F_g(i,j+2)+F_g(i,j)-2*F_g(i,j+1))/(h)^2;
        
        elseif (j == cols)
        %backward diff.
        d2Fdx2(i,j) = (F_g(i,j)+F_g(i,j-2)-2*F_g(i,j-1))/(h)^2;
        
        else
        %central difference
        d2Fdx2(i,j) = (F_g(i,j+1)+F_g(i,j-1)-2*F_g(i,j))/(h)^2;
        
        end
    end
end
grad_vec1 = d2Fdx2;
%% Compute the second part
[rows cols] = size(F_g);
for j = 1:cols
    for i = 1:rows
    %finite difference
        if (i == 1)
        %forward difference
        d2Fdt2(i,j) = (F_g(i+2,j)+F_g(i,j)-2*F_g(i+1,j))/(dt)^2;
        
        elseif (i == rows)
        %backward diff.
        d2Fdt2(i,j) = (F_g(i,j)+F_g(i-2,j)-2*F_g(i-1,j))/(dt)^2;
        
        else
        %central difference
        d2Fdt2(i,j) = (F_g(i+1,j)+F_g(i-1,j)-2*F_g(i,j))/(dt)^2;
        
        end
    end
end
grad_vec2 = d2Fdt2;

%% 
grad_vec = grad_vec1'+grad_vec2';

R_factor = R_factor_intensity * sqrt(dot(l_bottom,l_bottom))/...
    sqrt(dot(grad_vec,grad_vec)+10^(-50));

% R_factor_2 = R_factor_intensity * abs(max(l_bottom(:)))/...
%      (abs(max(grad_vec(:)))+10^(-50));
% grad2= -h*dt*l_bottom - R_factor_2 * h * dt * grad_vec; 
%%
grad = -h*dt*l_bottom - R_factor * h * dt * grad_vec;
