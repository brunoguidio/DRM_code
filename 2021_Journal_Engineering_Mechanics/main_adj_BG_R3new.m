clear%2D - Bruno - Summer 2019 
clear all; close all;  
tic
global sm; global h; global rho; global G; 
global T; global freq; global ne; global nn; global nb;
global time_intg_type; global dt; global beta; global gamma; 
global s_loc; global Ntimestep;global GK; global GM; global Keff; global GM_inv; 
global um; global Force_type; global Soil_profile; 
global BC; global GF_i; global h_best; global R_factor;global R_factor_intensity;

mkdir('graphs_F_g')
%mkdir('graphs_u')

%DEFINE PARAMETERS
Force_type = 'F5';            %choose: F1(smooth snake), F2(boomerang), F3(arc), F4(double arc), F5(Snake)
Soil_profile = 'Heterogeneous_6'; %choose: Homogeneous, Heterogeneous_3, Heterogeneous_6;
ns = 30;  %number os sensors

R_factor_intensity = 0.5;     %Ir
noise_level = 0; %[%]

Approach = 'OTD'; %OTD or DTO
%Set side 
sm = 60;      %side of mesh 
h_um = 0.5;         %side of square   
h_u = 1;

rho = 1000;    %mass density
if (strcmp(Soil_profile,'Homogeneous'))
    w_speed = [125]; 
elseif (strcmp(Soil_profile,'Heterogeneous_3'))
    %w_speed = [250;800;350];
    w_speed = [125;400;175];
elseif (strcmp(Soil_profile,'Heterogeneous_6'))
    %w_speed = [500;800;350;700;200;450];
    %w_speed = [250;400;175;350;100;225];
    w_speed = [400;450;300;350;200;250];
end
%G = w_speed.^2 * rho;      %Shear Modulus
for i = 1:length(w_speed)
    G(i,1) = w_speed(i,1)^2 * rho;
end
T = sm/max(w_speed)*4;

%Seismic_Input_Wave_type = 'Northridge';
Seismic_Input_Wave_type = 'Ricker';

freq = 10;  %frequency of Ricker pulse
h_best = 0.1;
Use_Conjugate_gradient_for_material = 'Y';

time_intg_type = 1; 
if (time_intg_type ==1)
    dt = 1e-3 ;%3e-04 %0.001;
    beta = 0.25; gamma = 0.5;   %Implicit time integration
    %mkdir('graphs_u_implicit')
else
    dt = 1e-5 ;%3e-04 %0.001;
    beta = 0; gamma = 0.5;      %explicit time integration
    %mkdir('graphs_u_explicit')
end

%% For u_m
h = h_um;
%Find elements and nodes
ne = (sm/h)^2;   %number of elements
nn = (sm/h+1)^2; %number of nodes

%Sensors location
nb = ((sm/h+1)-2)*2; %number of boundaries 

s = ceil((sm)/(ns+1));
for k = 1:ns
    if k == 1
        si = s -1;
        s_loc(k) = nn-nb-(sm/h)+((si/h)*k) ;
    else
        si = s;
        s_loc(k) = s_loc(k-1) + (si/h);
    end
end

[GK, GM, GF_i] = FEM_matrices( );
%%
if (strcmp(Seismic_Input_Wave_type,'Northridge'))
    Input_ugdot = load('Northridge_VT2_full.txt');           %read file, InputMotion: input motion
    Input_ugdot = Input_ugdot';    %convert columns to rows
    Input_ugdot = Input_ugdot(:);  %convert matrix columns to vector
    Input_ugdot = Input_ugdot/100; %cm/s => m/s
    truecp1 = Input_ugdot;
    time_interval_in_data = 0.005;
    truecp1 = fine_discretize(truecp1,time_interval_in_data,dt);
    [rows cols] = size(truecp1);
    num_cp = rows;
    
    for i = num_cp+1:num_cp*1.1
        truecp1(i) = 0.0;            %Add zero padding.
    end
    [rows cols] = size(truecp1);
    Ntimestep = rows;
    
    truecp = sparse(zeros(Ntimestep,sm/h+1));  %BG
    truecp(:,1) = truecp1;
    
    F_m = truecp_x_t(truecp);
    um = solve_u(F_m);

    %%%%%%
    F_g = sparse(zeros(Ntimestep,sm/h+1)); %BG
    A = 2e-3;
    for i = 1:Ntimestep
        t = i*dt;
        F_g(i,1) = ricker(freq,A,t);
    end

    F_g = force_x_t(F_g);
    F_g = max(truecp)/max(F_g)*F_g;
    
    u = solve_u(F_g);
%%%%%%%%% ------------------------------------ %%%%%%%%%%%%

elseif (strcmp(Seismic_Input_Wave_type,'Ricker'))

    Ntimestep = 1500; %T/dt;
    %Solve for um
    F_m = (zeros(Ntimestep,sm/h+1));
    A = 100;        %TRUE Ricker's amplitude
    for i = 1:Ntimestep
        t = i*dt;
        F_m(i,1) = ricker(freq,A,t);
    end
    
    F_m = force_x_t(F_m); 
    %F_m = flip(F_m,2); um = Function_add_noise_to_um(um,noise_level);
    
    um = solve_u(F_m);
    %um = Function_add_noise_to_um(um,noise_level);
    
end

%% For u
h = h_u;
%Find elements and nodes
ne = (sm/h)^2;   %number of elements
nn = (sm/h+1)^2; %number of nodes

%Sensors location
nb = ((sm/h+1)-2)*2; %number of boundaries 

s = ceil((sm)/(ns+1));
for k = 1:ns
    if k == 1
        si = s -1;
        s_loc(k) = nn-nb-(sm/h)+((si/h)*k) ;
    else
        si = s;
        s_loc(k) = s_loc(k-1) + (si/h);
    end
end


[GK, GM, GF_i] = FEM_matrices( );

%Solve for u - RICKER
    F_g = (zeros(Ntimestep,sm/h+1));  %%SPARSE
    A = 2e-3;
    for i = 1:Ntimestep
        t = i*dt;
        F_g(i,1) = ricker(freq,A,t);
    end

    F_g = force_x_t(F_g);
    %F_g = flip(F_g,2);
    u = solve_u(F_g);

%% Solving Adjoint Equation
if Approach == 'OTD'
l_bottom = solve_adj(u);

% Compute Gradient
grad = compute_grad(F_g,l_bottom);
%normalized gradient
grad_n = grad/abs(max(grad(:)));

% Updating the source profile
L = Function_errorf(um,u,F_g);
%L = Function_errorf_dto(um,u);

Src_D = -1 * grad;
Src_D_n = -1 * grad_n;
tol = 10^-18;
iteration_so_far = 1; %s
L_history(iteration_so_far) = L;
F_g_history{iteration_so_far,1} = F_g;
R_factor_history(iteration_so_far,1) = R_factor;
h_best_history(iteration_so_far,1) = h_best;
%
while (L > tol & iteration_so_far < 1000)
   iteration_so_far = iteration_so_far+1  
   F_g = F_g + h_best*Src_D_n';
   
   grad_prev = grad;
   Src_D_prev = Src_D; 
   
   u = solve_u(F_g);
   
   L = Function_errorf(um,u,F_g);
   %L = Function_errorf_dto(um,u);
   
   l_bottom = solve_adj(u);
   grad = compute_grad(F_g,l_bottom);

   %Search direction Src_D_Matr
   if (Use_Conjugate_gradient_for_material == 'Y')
     if (rem(iteration_so_far,5) == 0)
                Src_D = -1* grad;
     else
                Src_D = -1* grad + ...
                    dot(grad,grad)/dot(grad_prev,grad_prev)*...
                    Src_D_prev;
     end
   else
            Src_D = -1* grad;
   end
   % 
   
   Src_D_n = Src_D/abs(max(Src_D(:)));
   h_best = Newton_4_finding_h_best(F_g,Src_D_n);
   
   L_history(iteration_so_far) = L;
   F_g_history{iteration_so_far,1} = F_g;
   R_factor_history(iteration_so_far,1) = R_factor;
   h_best_history(iteration_so_far,1) = h_best;
   
    if (mod(iteration_so_far,10)==0)
    %clf;
    %Grid Matrix
    x = [1:sm/h+1];
    xm = [1:sm/h_um+1];
    y = [1:Ntimestep];
    ym= [1:Ntimestep];
    [xm,ym] = meshgrid(xm,ym);
   figure_F = figure(1); 
   subplot(1,2,1)
   hold on
   temp = [num2str(iteration_so_far),' -th iteration - TARGET']; 
   title(temp); 
   contourf(xm,ym,F_m,50,'LineColor','none')
   colorbar;
   caxis([-100 200])
   
   subplot(1,2,2)
   [x,y] = meshgrid(x,y);
   hold on
   temp = [num2str(iteration_so_far),' -th iteration - GUESS']; 
   title(temp); 
   contourf(x,y,F_g,50,'LineColor','none')
   colorbar;
   caxis([-100 200])
 
   filename =  ['graphs_F_g/F_g_at_' num2str(iteration_so_far) 'iteration.pdf']; 
            print(figure_F, '-r90', '-dpdf', filename);
            
   filename =  ['graphs_F_g/F_g_at_' num2str(iteration_so_far) 'iteration.fig'];  
            savefig(filename);
   
    %pause (1e-1)
    hold off; clf; close all;
    end
   
end

toc

elapsedTime = toc; 

filename = 'all_variables.mat';
save(filename)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
% Solving ADJ - using DTO
l_RBL = solve_adj_dto(u);

grad_dto = compute_grad_dto(F_g,l_RBL);
grad_dto_n = grad_dto/abs(max(grad_dto(:)));

% Updating the source profile
L = Function_errorf_dto(um,u);

Src_D = -1 * grad_dto;
Src_D_n = -1 * grad_dto_n;
tol = 10^-18;
iteration_so_far = 1; %s
L_history(iteration_so_far) = L;
F_g_history{iteration_so_far,1} = F_g;
R_factor_history(iteration_so_far,1) = R_factor;
h_best_history(iteration_so_far,1) = h_best;
%
while (L > tol & iteration_so_far < 1000)
   iteration_so_far = iteration_so_far+1  
   F_g = F_g + h_best*Src_D_n';
   
   grad_prev_dto = grad_dto;
   Src_D_prev = Src_D; 
   
   u = solve_u(F_g);
   
   L = Function_errorf_dto(um,u);
   
   l_RBL = solve_adj_dto(u);
   grad_dto = compute_grad_dto(F_g,l_RBL);
   
   %Search direction Src_D_Matr
   if (Use_Conjugate_gradient_for_material == 'Y')
     if (rem(iteration_so_far,5) == 0)
                Src_D = -1* grad_dto;
     else
                Src_D = -1* grad_dto + ...
                    dot(grad_dto,grad_dto)/dot(grad_prev_dto,grad_prev_dto)*...
                    Src_D_prev;
     end
   else
            Src_D = -1* grad_dto;
   end
   % 
   
   Src_D_n = Src_D/abs(max(Src_D(:)));
   h_best = Newton_4_finding_h_best_dto(F_g,Src_D_n);
   
   L_history(iteration_so_far) = L;
   F_g_history{iteration_so_far,1} = F_g;
   R_factor_history(iteration_so_far,1) = R_factor;
   h_best_history(iteration_so_far,1) = h_best;
   
    if (mod(iteration_so_far,50)==0)
    clf;
    %Grid Matrix
    x = [1:sm/h+1];
    xm = [1:sm/h_um+1];
    y = [1:Ntimestep];
    ym= [1:Ntimestep];
    [xm,ym] = meshgrid(xm,ym);
   figure_F = figure(1); 
   subplot(1,2,1)
   hold on
   temp = [num2str(iteration_so_far),' -th iteration - TARGET']; 
   title(temp); 
   contourf(xm,ym,F_m,50,'LineColor','none')
   colorbar;
   caxis([-100 200])
   
   subplot(1,2,2)
   [x,y] = meshgrid(x,y);
   hold on
   temp = [num2str(iteration_so_far),' -th iteration - GUESS']; 
   title(temp); 
   contourf(x,y,F_g,50,'LineColor','none')
   colorbar;
   caxis([-100 200])
 
   filename =  ['graphs_F_g/F_g_at_' num2str(iteration_so_far) 'iteration.pdf']; 
            print(figure_F, '-r90', '-dpdf', filename);
            
   filename =  ['graphs_F_g/F_g_at_' num2str(iteration_so_far) 'iteration.fig'];  
            savefig(filename);
   
   %pause (1e-1)
   hold off; clf; close all;
    end
   
end

toc

elapsedTime = toc; 

filename = 'all_variables_dto.mat';
save(filename)  
end