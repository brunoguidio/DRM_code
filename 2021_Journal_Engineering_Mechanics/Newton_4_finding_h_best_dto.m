function [ y ] = Newton_4_finding_h_best_dto (F_g,Src_D)

global um; global h_best;

epsilon_h = 10^-15; % preventing a zero denominator when computing an
%relative difference between h(i+1) and h(i) or
%preventing a zero denominator in gderiv.

RD_h = 1; %relative difference between the iterative solutions.
tol_h = 10^-5;

h_s(1) = h_best; %1/sqrt(dot(grad,grad)) %initial guess
dh =  (max(F_g(:))/max(Src_D(:)))*10^(-4); %dh for computing g derivative.

j = 1;

while (RD_h > tol_h & j <4)
    
    u = solve_u(F_g + Src_D' * (h_s(j)+dh));       
    %Misfit1 = Function_errorf(um,u,F_g);
    Misfit1 = Function_errorf_dto(um,u);

    u = solve_u(F_g + Src_D' * (h_s(j)-dh));       
    %Misfit2 = Function_errorf(um,u,F_g);
    Misfit2 = Function_errorf_dto(um,u);

    u = solve_u(F_g + Src_D' * (h_s(j)));       
    %Misfit3 = Function_errorf(um,u,F_g);
    Misfit3 = Function_errorf_dto(um,u);
    
    h_deriv_mat = (Misfit1-Misfit2)/(2*dh);
    
    h_double_deriv_mat = (Misfit1+Misfit2-2*Misfit3)/(dh)^2;
    
    h_s(j+1) = h_s(j) -h_deriv_mat/(h_double_deriv_mat+epsilon_h);
    
    if (isnan(h_s(j+1)))
        h_s(j+1) = h_s(j)/2;
    end      
        
    if (h_s(j+1)<0)
        h_s(j+1) = -1*h_s(j+1); 
    end 
    
    RD_h = abs(h_s(j+1) - h_s(j))...
        /abs(h_s(j)+epsilon_h)  *100;
    RD_h_arr(j) = RD_h;
    j = j+1;
end

y = h_s(j); %*0.9;
    

