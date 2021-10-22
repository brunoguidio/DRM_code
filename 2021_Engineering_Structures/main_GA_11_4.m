%Example 11 - test 4
warning('off','all')
%5 parameters - 2 source and 3 structural
tic
clear all; close all;
global P_true; global freq_true; global mv; global x0;
global E_true_1;global E_true_2;global E_true_3;
%global E_true_4;global E_true_5;
%global E_true_6;global E_true_7;global E_true_8;global E_true_9;

global ne; global ns; global u_true; global iteration_so_far; global L
global FGM; global dt; %global FGK; global Keff; global Keff_inv;  

global error_mat_each_it; 
global GN; global PS;

global P_population_history; global freq_population_history; 
global E1_population_history;global E2_population_history;
global E3_population_history;%global E4_population_history;
%global E5_population_history;global E6_population_history;
%global E7_population_history;global E8_population_history;
%global E9_population_history;%global E10_population_history;

global N_GS_Load; 
N_GS_Load = 1; 

mv = 20; x0 = 50;
 
P_true = 100; freq_true = 20;

E_true_1 = 1.8;%*10^10;
E_true_2 = 2.5;%*10^10;
E_true_3 = 2.5;%*10^10;
% E_true_4 = 2.5;%*10^10;
% E_true_5 = 2.5;%*10^10;
% E_true_6 = 1.8;%*10^10; %***
% E_true_7 = 2.5;%*10^10;
% E_true_8 = 2.5;%*10^10; 
% E_true_9 = 2.5;%*10^10;
%E_true_10 = 2.5*10^10;

GN = 50; 
PS = 100; 
ns = 45; 

L = 100;
ne = 180;
 
error_mat_each_it = zeros(GN,PS);

P_population_history = zeros(N_GS_Load,GN,PS);
freq_population_history = zeros(N_GS_Load,GN,PS);
E1_population_history = zeros(N_GS_Load,GN,PS);
E2_population_history = zeros(N_GS_Load,GN,PS);
E3_population_history = zeros(N_GS_Load,GN,PS);
% E4_population_history = zeros(N_GS_Load,GN,PS);
% E5_population_history = zeros(N_GS_Load,GN,PS);
% E6_population_history = zeros(N_GS_Load,GN,PS);
% E7_population_history = zeros(N_GS_Load,GN,PS);
% E8_population_history = zeros(N_GS_Load,GN,PS);
% E9_population_history = zeros(N_GS_Load,GN,PS);
%E10_population_history = zeros(N_GS_Load,GN,PS);

iteration_so_far = 0; 
FGM = global_mass_matrix (ne);
dt = 0.001; 
%FGM = sparse(FGM);

u_true = CB_ML_5 (P_true,freq_true,...
                  E_true_1,E_true_2,E_true_3,ne,ns);
%%              
fprintf('starting error func \n');
fitfunc = @error_movingload;

dev_source = 0.5; 
dev_struct = 0.05;

low_P = P_true*(1-dev_source);   up_P = P_true*(1+dev_source);
low_freq = freq_true*(1-dev_source);  up_freq = freq_true*(1+dev_source);
LP = zeros(1,N_GS_Load); LP(:) = low_P ;
LF = zeros(1,N_GS_Load); LF(:) = low_freq ; 
UP = zeros(1,N_GS_Load); UP(:) = up_P ;
UF = zeros(1,N_GS_Load); UF(:) = up_freq ; 

E_min = min([E_true_1,E_true_2,E_true_3]);
E_max = max([E_true_1,E_true_2,E_true_3]);
LE = E_min*(1-dev_struct); UE = E_max*(1+dev_struct);
%%
opts = gaoptimset('Generations',GN, 'PopulationSize', PS);
[soln fval] = ga(fitfunc,N_GS_Load*2+3,[],[],[],[],...
                 [LP LF LE LE LE],...
                 [UP UF UE UE UE],[],[],opts);
             
             
             
             
             
             
soln 
toc

elapsedTime = toc; 

filename = 'all_variables.mat';
save(filename)