function u = CB_ML_5 (P, freq, E_1,E_2,E_3, n, n_sensors)
%%
P = P_true;freq = freq_true;
E_1 = E_true_1;E_2 = E_true_2;E_3 = E_true_3;
n = ne;n_sensors = ns;
% % 
% mkdir('graphs_force')
% 
% P = soln(1);freq = soln(2);
% E_1 = soln(3);E_2 = soln(4);E_3 = soln(5);
% n = ne;n_sensors = ns;

Ntimestep = 1000;
global L;
h = L/n;
global mv;
global x0;
global dt; 
global FGM; 
% global FGK; % global Keff; % global Keff_inv; 

FGK = global_stiffness_matrix (n,E_1,E_2,E_3); 
Keff = FGM + FGK*dt*dt*0.25; 
Keff_inv = inv(Keff); 
Keff_inv =sparse(Keff_inv);
FGK = sparse(FGK);

%fprintf('starting CB_ML \n');
t = 0;

GF = sparse(zeros (5*n+2,1));

%GLOBAL Force Matrix
for nd = 1:n
    
    x1 = (nd-1) * h;
    x2 = x1 + h;
    
    
    for i = 1:4
        
        F3_LE_cmp = @(xi)(sh_3_func(i,xi))*(h/2).*mult_mov_load_guess(xi, x1, x2, mv, P, freq, x0, t);
          
        F3_LE(i,1) = quad(F3_LE_cmp,-1,1);
        
    end
    
    GF (3*nd-2:3*nd+1,1)= GF(3*nd-2:3*nd+1,1)+ F3_LE;
%     y (3*nd-2:3*nd+1,1) = [mult_mov_load_guess(-1,x1,x2,mv,P,freq,x0,t);...
%                                 mult_mov_load_guess(-1/3,x1,x2,mv,P,freq,x0,t);...
%                                 mult_mov_load_guess(1/3,x1,x2,mv,P,freq,x0,t);...
%                                 mult_mov_load_guess(1,x1,x2,mv,P,freq,x0,t)];
                                                   
end

%GF(1:3*n+1,1) = GF(1:3*n+1,1).*y; 
GF(3*n+1,:)=[]; GF(2*n+1,:)=[];
GF(n+1,:)=[];  GF(1,:)=[];

u_curt = 0;
udot_curt = 0;
uddot_curt = FGM\GF;

    spacing = round((((3*n+1)-n_sensors)/(n_sensors+1))+1);
    
% sensors_location(1) = (spacing-1) * L/(3*n);
% for tn = 2: n_sensors
% sensors_location(tn) = sensors_location(tn-1) + (spacing) * L/(3*n);             
% end
%fprintf('loop \n');
%tic
%%
for ts = 1: Ntimestep
    
    t = ts*dt;
    u_prev = u_curt;
    udot_prev = udot_curt;
    uddot_prev = uddot_curt;
    
    %Build the RHS load vector
    GF = sparse(zeros (5*n+2,1));
    
    %GLOBAL Force Matrix
    for nd = 1:n
        
        x1 = (nd-1) * h;
        x2 = x1 + h;
        
        for i = 1:4
           
            F3_LE_cmp = @(xi)(sh_3_func(i,xi))*(h/2).*mult_mov_load_guess(xi, x1, x2, mv,P,freq,x0, t);
            
            F3_LE(i,1) = quad(F3_LE_cmp,-1,1);
        
        end
        
        GF (3*nd-2:3*nd+1,1)= GF(3*nd-2:3*nd+1,1)+ F3_LE;
%         y (3*nd-2:3*nd+1,1) = [mult_mov_load_guess(-1,x1,x2,mv,P,freq,x0,t);...
%                                 mult_mov_load_guess(-1/3,x1,x2,mv,P,freq,x0,t);...
%                                 mult_mov_load_guess(1/3,x1,x2,mv,P,freq,x0,t);...
%                                 mult_mov_load_guess(1,x1,x2,mv,P,freq,x0,t)];

                             
    end
 %   GF(1:3*n+1,1) = GF(1:3*n+1,1).*y;
    GF(3*n+1,:)=[]; GF(2*n+1,:)=[];
    GF(n+1,:)=[];  GF(1,:)=[];
    
    RHS = GF - FGK*(u_prev + udot_prev*dt + 0.25*uddot_prev*dt^2);
    
    %Solve for acceleration at the current (i)-th time step
    %uddot_curt = Keff\RHS;
    uddot_curt = Keff_inv * RHS; %CJ 2016/7/21
    
    %Updating u and udot at the current (i)-th time step
    u_curt = u_prev + udot_prev*dt + 0.5*0.5*(uddot_prev + uddot_curt)*dt*dt;
    udot_curt = udot_prev + 0.5*(uddot_prev + uddot_curt)*dt;
    
    u_total = [0;u_curt(1:n-1);0;u_curt(n:2*n-2);0;u_curt(2*n-1:3*n-3);0];
    
    %CJ
    %u_bruno(ts,:) = [0;u_curt(1:n-1);0;u_curt(n:2*n-2);0;u_curt(2*n-1:3*n-3);0];
    
    %this u is for sensors_arr
    u(ts,:) = u_total(spacing:spacing:3*n+1);
    
%     if (mod(ts,10)==0)
%       
%         figure(1); clf;
%         hold on
%         temp = [num2str(t),' second'];
%         title(temp);
% 
%         x_arr = linspace(0,L,3*n+1);
%         plot(x_arr,u_total)
%         %ylim([-2e-5,2e-5])
%         pause(1e-1)
%     end

%     if (mod(ts,20)==0)
%         fig2 = figure(2); clf;
%         hold on
%         x_arr = linspace(0,L,3*n+1);
%         [hAx,hLine1,hLine2] = plotyy(x_arr,u_total,x_arr,y);
%         grid on;
%         axis(hAx(1),[0,L,-2e-4,2e-4])
%         axis(hAx(2),[0,L,-250,250])
%         set(hAx(1),'Ylim',[-2e-4,2e-4])
%         xlabel('x [m]','Interpreter','Latex','FontSize',12);
%         ylabel(hAx(1),'wave response [m]','Interpreter','Latex','FontSize',12);
%         ylabel(hAx(2),'Source Function [N/m]','Interpreter','Latex','FontSize',12);
%         legend('Wave responses','Source Function','Location','SouthEast','FontSize',10);
% 
%         pause (1e-10)
%         
%     end


%         if (mod(ts,10)==0)
%            fig3 = figure(3); clf;
%            hold on; 
%            x_arr = linspace(0,L,3*n+1);
%            yyaxis left
%            plot(x_arr,u_total); grid on
%            ylim([-1.5e-4 1.5e-4])
%            xlabel('$x$ [m]','Interpreter','Latex','FontSize',12);
%            ylabel('wave response [m]','Interpreter','Latex','FontSize',12);
%            
%            yyaxis right
%            plot(x_arr,y,'--')
%            ylim([-200 200])
%            
%            ylabel('Source Function [N/m]','Interpreter','Latex','FontSize',12);
%            legend('$w$','$q$','Interpreter','Latex','Location','NorthEast','FontSize',10);
%            
%            filename =  ['graphs_force/Wave_Force_at_' num2str(t) 'seconds.pdf'];   
% %     
%            print(fig3, '-r90', '-dpdf', filename); 
%         end


end
%toc

