P = P_true;freq = freq_true;
E_1 = E_true_1;E_2 = E_true_2;E_3 = E_true_3;
n = ne;n_sensors = ns;

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
        
        F3_LE_cmp = @(xi)(sh_3_func(i,xi))*(h/2);%.*mult_mov_load_guess(xi, x1, x2, mv, P, freq, x0, t);
          
        F3_LE(i,1) = quad(F3_LE_cmp,-1,1);
        
    end
    
    GF (3*nd-2:3*nd+1,1)= GF(3*nd-2:3*nd+1,1)+ F3_LE;
    y (3*nd-2:3*nd+1,1) = [mult_mov_load_guess(-1,x1,x2,mv,P,freq,x0,t);...
                                mult_mov_load_guess(-1/3,x1,x2,mv,P,freq,x0,t);...
                                mult_mov_load_guess(1/3,x1,x2,mv,P,freq,x0,t);...
                                mult_mov_load_guess(1,x1,x2,mv,P,freq,x0,t)];
                            
                            
    ys (3*nd-2:3*nd+1,1) = [load_q_mov_guess(-1,x1,x2,mv,P,freq,x0,t);...
                                load_q_mov_guess(-1/3,x1,x2,mv,P,freq,x0,t);...
                                load_q_mov_guess(1/3,x1,x2,mv,P,freq,x0,t);...
                                load_q_mov_guess(1,x1,x2,mv,P,freq,x0,t)];
                                                   
end

GF(1:3*n+1,1) = GF(1:3*n+1,1).*y; 
GF(3*n+1,:)=[]; GF(2*n+1,:)=[];
GF(n+1,:)=[];  GF(1,:)=[];

u_curt = 0;
udot_curt = 0;
uddot_curt = FGM\GF;

    spacing = round((((3*n+1)-n_sensors)/(n_sensors+1))+1);
    
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
           
            F3_LE_cmp = @(xi)(sh_3_func(i,xi))*(h/2);%.*mult_mov_load_guess(xi, x1, x2, mv,P,freq,x0, t);
            
            F3_LE(i,1) = quad(F3_LE_cmp,-1,1);
        
        end
        
        GF (3*nd-2:3*nd+1,1)= GF(3*nd-2:3*nd+1,1)+ F3_LE;
        y (3*nd-2:3*nd+1,1) = [mult_mov_load_guess(-1,x1,x2,mv,P,freq,x0,t);...
                                mult_mov_load_guess(-1/3,x1,x2,mv,P,freq,x0,t);...
                                mult_mov_load_guess(1/3,x1,x2,mv,P,freq,x0,t);...
                                mult_mov_load_guess(1,x1,x2,mv,P,freq,x0,t)];
                            
        ys (3*nd-2:3*nd+1,ts+1) = [load_q_mov_guess(-1,x1,x2,mv,P,freq,x0,t);...
                                load_q_mov_guess(-1/3,x1,x2,mv,P,freq,x0,t);...
                                load_q_mov_guess(1/3,x1,x2,mv,P,freq,x0,t);...
                                load_q_mov_guess(1,x1,x2,mv,P,freq,x0,t)];
                                                   

                             
    end
    GF(1:3*n+1,1) = GF(1:3*n+1,1).*y;
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
end

%%
        fig = figure(1); clf;
        hold on; grid on;  
%         temp = [num2str(t),' second'];
%         title(temp);
        x_arr = linspace(0,L,3*n+1);
        plot(x_arr,ys(:,1))
        plot(x_arr,ys(:,end),'--')
        
        xlabel('$x$ [m]','Interpreter','Latex','FontSize',12);
        ylabel('$H(x,t)$','Interpreter','Latex','FontSize',12);
        
        legend('$t$ = 0 s','$t$ = 1 s','Location','NorthEast','FontSize',10,'Interpreter','latex');
        set(gcf, 'Position', [200,200,600,300])
        
        xlim([40 80])
        %ylim([0 3])
        
        filename =  ['H_spatial_variation.pdf'];   
%     
        print(fig, '-r90', '-dpdf', filename); 