clear all; clc; 
filename = 'u_true.mat';
load(filename);

u_true_40 = u_true(:,18);
u_true_60 = u_true(:,28);
t = [0.001:0.001:1];

filename = 'u_soln.mat';
load(filename);

u_soln_40 = u(:,18);
u_soln_60 = u(:,28);
%%
fig_handle = figure (1);
hold on
plot(t,u_true_40,'b')
plot(t,u_soln_40,'r--','LineWidth',1.5)
legend('$w^m$','$w$','Interpreter','Latex','Location','SouthEast','FontSize',10);
        
    ylabel('wave response [m]','Interpreter','Latex',...
        'FontSize',12);
    xlabel('time [s]','Interpreter','Latex',...
        'FontSize',12);
       filename =  ['wave_response_40.pdf'];   
%     
       print(fig_handle, '-r90', '-dpdf', filename); 
%%
fig_handle = figure(2);
hold on
plot(t,u_true_60,'b')
plot(t,u_soln_60,'r--','LineWidth',1.5)
legend('$w^m$','$w$','Interpreter','Latex','Location','SouthEast','FontSize',10);
        
    ylabel('wave response [m]','Interpreter','Latex',...
        'FontSize',12);
    xlabel('time [s]','Interpreter','Latex',...
        'FontSize',12);
       filename =  ['wave_response_60.pdf'];   
%     
       print(fig_handle, '-r90', '-dpdf', filename); 

%%
filename = 'u_total_true.mat';
load(filename);
u_total_true = u_total;

filename = 'u_total_soln.mat';
load(filename);
u_total_soln = u_total;

x_arr = linspace(0,100,541);
fig_handle = figure(3);
hold on
plot(x_arr,u_total_true,'b')
plot(x_arr,u_total_soln,'r--','LineWidth',1.5)
legend('u^m','u','Location','SouthEast','FontSize',10);
        
    ylabel('wave response [m]','Interpreter','Latex',...
        'FontSize',12);
    xlabel('$x$ [m]','Interpreter','Latex',...
        'FontSize',12);
    
       filename =  ['wave_response_1second.pdf'];   
%     
       print(fig_handle, '-r90', '-dpdf', filename); 
       
       %%
U_B = fft(u_true_40); 
U_B = abs(U_B);

U_F = fft(u_true_60); 
U_F = abs(U_F);