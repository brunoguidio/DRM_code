function [L] = Function_errorf(um,u,F_g)
global h; global dt;global R_factor;

[row col] = size(um); 
    z = zeros(row,col);
    for k = 1:row
        for j = 1:col
            z(k,j) = (um(k,j)-u(k,j))^2;
        end 
    end 
   L = sum(sum(z))*dt;
      
%% Compute the Regularization

Reg_term = 0;
[rows cols] = size(F_g);
for i = 1:rows
    for j = 1:cols
        if (j == 1)
        %forward difference
        dFdx = (F_g(i,j+1)-F_g(i,j))/(h);   
        
        elseif (j == cols)
        %backward diff.
        dFdx = (F_g(i,j)-F_g(i,j-1))/(h);
        
        else
        %central difference
        dFdx = (F_g(i,j+1)-F_g(i,j-1))/(2*h);
        end
        Reg_term = Reg_term+(dFdx)^2 * h * dt; %TN
    end
end
%% Compute second part 
[rows cols] = size(F_g);
for j = 1:cols
    for i = 1:rows
        if (i == 1)
        %forward difference
        dFdt = (F_g(i+1,j)-F_g(i,j))/(dt);   
        
        elseif (i == rows)
        %backward diff.
        dFdt = (F_g(i,j)-F_g(i-1,j))/(dt);
        
        else
        %central difference
        dFdt = (F_g(i+1,j)-F_g(i-1,j))/(2*dt);
        end
        Reg_term = Reg_term+(dFdt)^2 * h * dt; %TN
    end
end


%%
Reg_term = 0.5* R_factor*Reg_term;
%y = y+Reg_term; %with regularization
L = L+Reg_term;