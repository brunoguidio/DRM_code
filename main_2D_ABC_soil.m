% 2D Plane-Strain and DRM  % Written by Bruno Guidio

clear all; clc;
global L; global H; global h; global NV; global NH; 
global rho; global G; global E; global nu; 
global Vs; global Vp; global Mod_p; global Zs; global Zp;  
global N_Ele; global N_Node; global N_DOF;
global dt; global Ntimestep; 
global GK; global GM; global GC; global Keff; global Meff;
global Soil_profile; global x_arr; global y_arr;
global s_loc; global um; global h_best; 
global corner_of_b; global size_of_b; global corner_of_i; global size_of_i;
global corner_of_e; global size_of_e; global number_of_layers; global number_of_inclusions
global Li; global damping_a; global small_domain;
global n_i; global n_b; global n_e; global n_ge; global n_e_out;  global n_ge1; global n_ge2;
global L1; global GM_t; global GC_t; global Keff_t;

%mkdir('graphs_u'); mkdir('graphs_u2')
% mkdir('graphs_ux')
% mkdir('graphs_uy')
load('MyColormap.mat')

% Dimensions of a cantilever beam. 
L = 80;                             %Length 
H = 40;                             %Height

h = 1;                              %Element size 

size_of_small_domain = [40,20];     %[length,height]
number_of_layers = 4;
number_of_inclusions = 2;

dt = 0.001;
Ntimestep = 1500;
Ricker_freq = 2;                   %Hz

NV = H/h;                           %elements in the vertical axis
NH = L/h;                           %elements in the horizontal axis
N_Ele = NV * NH;                    %number of elements
N_Node = (2*NV+1)*(2*NH+1);         %number of nodes
N_DOF = 2*N_Node;                   %number of degrees of freedom before tailoring.

%Loading information 
%Point loading at the right end at the center height. 
Point_loading_ct_right_x = 0; %-100% -100; %N/m;
Point_loading_ct_right_y = 0; %; -100; %N/m;

Point_loading_x = point_loading(Ricker_freq);

% For postprocessing. 
x_arr = [0:1:2*(NH)];
x_arr = x_arr *h/2;

y_arr = [0:1:2*(NV)];
y_arr = y_arr *h/2;

%% Material property

%Soil
rho = 1500;   %1800!!                               % mass density - kg/m^3

% Vs = [200;200;200;200;500;800];                   % Shear wave speed [m/s]
% Vp = [380;380;380;380;980;1600];

Vs = [200;150;150;100;500;800]; %[200;200;200;200;500];                       % Shear wave speed [m/s]
Vp = [400;300;300;200;1000;1600]; %[380;380;380;380;850]; 

nu = (Vp.^2-2*Vs.^2)./(2*(Vp.^2-Vs.^2));
G = rho.*(Vs.^2);                                   % Shear modulus [N/m2]
E = (1 + nu).*(2*G);                                % Young's modulus [N/m2]
Mod_p = 2*G.*(1 - nu)./(1 - 2*nu);                  % P-wave modulus 
Zp = rho.*Vp;
Zs = rho.*Vs;   

%% Inclusions
Li{1} = [40,50; 24,28];
Li{2} = [28,33; 27,40];
% Li{1} = [0.5*L,0.6*L;0.68*H,0.76*H];
% Li{2} = [0.38*L,0.43*L;0.74*H,H];

%% Step 1: Enter connectivity array.
[Connectivity] = connectivity_matrix ();

%% Step 2 Building Matrices
Count_Inclusions = 'N';
[GK,GM,GC,GF] = building_matrices_hetero (Connectivity,Count_Inclusions);

%% Step 3 Body force

GF_m = zeros(N_DOF,Ntimestep+1);
point_load_loc = [L/2,5];         %body arc

GF_m(N_Node+(2*NH+1)*(2*point_load_loc(2)/h)+1+(2*point_load_loc(1)/h),:) = ...
                      GF(N_Node+(2*NH+1)*(2*point_load_loc(2)/h)+1+(2*point_load_loc(1)/h),:) + Point_loading_x';
%% Step 4 Effective force
size_of_e = size_of_small_domain;                       %[width,height]
corner_of_e = [(L-size_of_e(1))/2,H-size_of_e(2)];      %[x,y]

corner_of_b = [corner_of_e(1)+2,corner_of_e(2)+2];      %[x,y]             % ATTENTION distance to gamma_b
size_of_b = [size_of_e(1)-4,size_of_e(2)-2];           %[width,height]

corner_of_i = [corner_of_b(1)+h/2,corner_of_b(2)+h/2];  %[x,y]
size_of_i = [size_of_b(1)-h,size_of_b(2)-h/2];          %[width,height]

F_eff_history2 = solve_F_eff_2D_bigger(GF_m);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Bigger Domain with inclusions

% Building Matrices with inclusions
Count_Inclusions = 'Y';
[GK,GM,GC,GF] = building_matrices_hetero (Connectivity,'Y');
%% Step 5 : Solve u using bigger domain and GF_m

[u_curt_history_big] = solve_u_DRM (GF_m);
%um = u_curt_history_big(s_loc,:);

nodes_u = [];
for dep = corner_of_e(2):h/2:H
    
row = [2*corner_of_e(1)+1:2*(corner_of_e(1)+size_of_e(1))+1] + ((2*L/h+1)*(2*dep/h)); 

nodes_u = [nodes_u,row];
end
nodes_u = [nodes_u,(nodes_u+N_Node)];

u_curt_history_2 = u_curt_history_big(nodes_u,:);

%% Small Domain with inclusions

L = size_of_small_domain(1);  %400                       %Length 
H = size_of_small_domain(2);  %120                       %Height
number_of_layers = 2;
        
NV = H/h;                           %nodes in the vertical axis
NH = L/h;                           %nodes in the horizontal axis
N_Ele = NV * NH;                    %number of elements
N_Node = (2*NV+1)*(2*NH+1);         %number of nodes
N_DOF = 2*N_Node;                   %number of degrees of freedom before tailoring.

% For postprocessing. 
x_arr = [0:1:2*(NH)];
x_arr = x_arr *h/2;

y_arr = [0:1:2*(NV)];
y_arr = y_arr *h/2;

%% Material property

%Soil
Vs = Vs(3:6);
Vp = Vp(3:6);
nu = (Vp.^2-2*Vs.^2)./(2*(Vp.^2-Vs.^2));
G = rho.*(Vs.^2);                   % Shear modulus [N/m2]
E = (1 + nu).*(2*G);                % Young's modulus [N/m2]
Mod_p = 2*G.*(1 - nu)./(1 - 2*nu);  % P-wave modulus 
%Vp = sqrt(Mod_p./rho);              % P-wave wave speed [m/s]
Zp = rho.*Vp;
Zs = rho.*Vs; 

%% Step 1: Enter connectivity array.
[Connectivity] = connectivity_matrix ();

%% Step 2 Building Matrices
small_domain = 'Y';
%damping_a = 30000;
%damping_a = [25000,25000,25000];
%damping_a = [150000,150000,150000];

Li{1} = [Li{1}(1,1)-20,Li{1}(1,2)-20;Li{1}(2,1)-20,Li{1}(2,2)-20];
Li{2} = [Li{2}(1,1)-20,Li{2}(1,2)-20;Li{2}(2,1)-20,Li{2}(2,2)-20];
Count_Inclusions = 'Y';
[GK,GM,GC,GF] = building_matrices_hetero (Connectivity,Count_Inclusions);

%% Step 3 Find n_i, n_b, n_ge, and n_e again
%corner_of_b = [corner_of_b(1)-10,corner_of_b(2)-5];           %[x,y]
corner_of_b = [corner_of_b(1)-corner_of_e(1),corner_of_b(2)-corner_of_e(2)]; 
%size_of_b = [size_of_e(1)-20,size_of_e(2)-5];           %[width,height]

%corner_of_i = [corner_of_i(1)-10,corner_of_i(2)-5];          %[x,y]
corner_of_i = [corner_of_b(1)+h/2,corner_of_b(2)+h/2];
%size_of_i = [size_of_b(1)-2,size_of_b(2)-1];                %[width,height]
nodes_small_domain();

spac_meter = 1;  
spac = 2*spac_meter;                % 2 means two nodes, 1 m
%Sensors location
s_loc = [2*corner_of_b(1):spac:2*(corner_of_b(1)+size_of_b(1))] + ((2*L/h+1)*(2*H/h)+1); 
s_loc = [s_loc,(s_loc+N_Node)];
%% Step 5 : Solve u using bigger domain and GF_m
[u_curt_history_small] = solve_u_DRM (F_eff_history2);
%um_DRM = u_curt_history_small(s_loc,:);

%
error_ui = sum(abs(u_curt_history_small(n_i,:)-u_curt_history_2(n_i,:)).^2)/sum(abs(u_curt_history_2(n_i,:)).^2)*100
error_ui2 = sum(abs(u_curt_history_small(n_i,:)-u_curt_history_2(n_i,:)).^1)/sum(abs(u_curt_history_2(n_i,:)).^1)*100

%%
%um_abc = u_curt_history(s_loc,:);
load('u_curt_history_PML.mat')
um = u_curt_history_PML(s_loc,:);

%um = u_curt_history_2(s_loc,:);
%%

%F_eff_history_g = sparse(zeros(size(F_eff_history2)));
F_eff_history_g = sparse(zeros(2*N_Node,Ntimestep+1));

[u_curt_history_g] = solve_u_DRM (F_eff_history_g);
u = u_curt_history_g(s_loc,:);

%% Objective functional

L_error = Function_errorf_dto(um,u);

%% Solving Adjoint Equation
Keff = sparse((4/dt^2)*GM+(2/dt)*GC+GK);
Keff_t = Keff';
%Keff_t_inv = inv(Keff_t);
GM_t = GM';
%GM_t_inv = GM_inv';
GC_t = GC';
%tic
l_eff = solve_adj_DRM(u);
%%
grad = (l_eff);

grad_n = grad/max(abs(grad(:)));

%% Update Search direction

iteration_so_far = 1;

Src_D = -1 * grad;
Src_D_n = -1 * grad_n;
tol = L_error*10^-7;
L_error_history(iteration_so_far) = L_error;
h_best = 0.1; Use_Conjugate_gradient_for_material = 'Y';
h_best_history(iteration_so_far,1) = h_best;


%%
u_curt_history = u_curt_history_PML;
%u_curt_history = u_curt_history_2;
error_ui = sum(abs(u_curt_history_g(n_i,:)-u_curt_history(n_i,:)).^2)/sum(abs(u_curt_history(n_i,:)).^2)*100;
error_ui2 = sum(abs(u_curt_history_g(n_i,:)-u_curt_history(n_i,:)).^1)/sum(abs(u_curt_history(n_i,:)).^1)*100;
error_ub = sum(abs(u_curt_history_g([n_b],:)-u_curt_history([n_b],:)).^2)/sum(abs(u_curt_history([n_b],:)).^2)*100;
error_uge = sum(abs(u_curt_history_g(n_ge,:)-u_curt_history(n_ge,:)).^2)/sum(abs(u_curt_history(n_ge,:)).^2)*100;
        
error_ui_history(iteration_so_far) = error_ui;
error_ui2_history(iteration_so_far) = error_ui2;
error_ub_history(iteration_so_far) = error_ub;
error_uge_history(iteration_so_far) = error_uge;

%%
while (L > tol & iteration_so_far < 5)
   
   iteration_so_far = iteration_so_far+1
   F_eff_history_g([n_b,n_ge],:) = F_eff_history_g([n_b,n_ge],:) + h_best*Src_D_n;
   
   [u_curt_history_g] = solve_u_DRM (F_eff_history_g);
   u = u_curt_history_g(s_loc,:);
   
   % <<<<<>>>>>> %
   error_ui = sum(abs(u_curt_history_g(n_i,:)-u_curt_history(n_i,:)).^2)/sum(abs(u_curt_history(n_i,:)).^2)*100
   error_ui2 = sum(abs(u_curt_history_g(n_i,:)-u_curt_history(n_i,:)).^1)/sum(abs(u_curt_history(n_i,:)).^1)*100
   % <<<<<>>>>>> % 
   
   L_error = Function_errorf_dto(um,u);
   
   grad_prev = grad;
   Src_D_prev = Src_D; 
       
   l_eff = solve_adj_DRM(u);
   grad = l_eff;

%   Src_D = -1* grad;
    %Search direction Src_D_Matr
    if (Use_Conjugate_gradient_for_material == 'Y')
        if (rem(iteration_so_far,20) == 0)
            Src_D = -1* grad;
        else
            Src_D = -1* grad + dot(grad,grad)/dot(grad_prev,grad_prev)*Src_D_prev;
        end
    else
            Src_D = -1* grad;
    end
   
   
   Src_D_n = Src_D/max(abs(Src_D(:)));
   
  
   h_best = Newton_4_finding_h_best(F_eff_history_g,Src_D_n);
  
   L_error_history(iteration_so_far) = L_error;
   h_best_history(iteration_so_far,1) = h_best;
  
    error_ub = sum(abs(u_curt_history_g([n_b],:)-u_curt_history([n_b],:)).^2)/sum(abs(u_curt_history([n_b],:)).^2)*100;
    error_uge = sum(abs(u_curt_history_g(n_ge,:)-u_curt_history(n_ge,:)).^2)/sum(abs(u_curt_history(n_ge,:)).^2)*100;
        
    error_ui_history(iteration_so_far) = error_ui;
    error_ui2_history(iteration_so_far) = error_ui2;
    error_ub_history(iteration_so_far) = error_ub;
    error_uge_history(iteration_so_far) = error_uge;  
   
    if (mod(iteration_so_far,50)==0)
      filename = ['all_variables_' num2str(iteration_so_far) '.mat'];
      save(filename)  
    end
   
end
