function [u_sensors] = solve_u (F_matrix)

global sm; global h; global rho; global G;
global T; global freq; global ne; global nn; global nb;
global time_intg_type; global dt; global beta; global gamma; 
global s_loc; global Ntimestep;
global GK; global GM; global Keff; global GM_inv; global MM; global BC
global GF_i


F_vector = sparse(zeros(Ntimestep,nn-nb));
F_vector(1:Ntimestep,1:sm/h+1) = F_matrix;



%GF= -1*GF_i*ricker(freq,A,t);   % 
GF= -1*GF_i*0;

%%%%% Boundary Conditions %%%%%%%%%
%GF = GF(BC,1);

%Newmark Time Integration
Keff = GM + GK*dt*dt*0.25;

GF = sparse(GF);
%Keff = sparse(Keff);

u_curt = zeros(length(GF),1);
udot_curt = zeros(length(GF),1);
uddot_curt = GM_inv * GF;
%%
for ts = 1: Ntimestep
    
    t = ts*dt; 
    u_prev = u_curt;
    udot_prev = udot_curt;
    uddot_prev = uddot_curt;
    
    %Build the RHS load vector
     %Global Force Vector
        
    GF= -1*GF_i.*F_vector(ts,:)';   % STRESS??? *

    if (time_intg_type ==1)
        %Implicit time integration
    RHS = GF - GK*(u_prev + udot_prev*dt + 0.25*uddot_prev*dt^2);
    %Solve for acceleration at the current (i)-th time step
    uddot_curt = Keff\RHS;
    else       
        %explicit time integration.
        RHS = GF - GK * (u_prev);
        uddot_curt = GM_inv * RHS;
    end
    
    %Updating u and udot at the current (i)-th time step
    u_curt = u_prev + udot_prev * dt + ...
        0.5 * (uddot_prev * (1.0-2.0*beta) + uddot_curt* 2.0*beta)*dt^2;
    
    udot_curt = udot_prev + ...
        (uddot_prev * (1.0-gamma) + uddot_curt* gamma)*dt;
    
    %u_curt_history(:,ts) = u_curt;
    
    for i = 1:length(s_loc)
        u_sensors (i,ts) = u_curt(s_loc(i),1);
    end
    
%     %Convert solution u_curt to matrix
%     W = zeros(length(u_curt),1);
%     W_prev = u_curt;
%     b=[0];
%     for i =2:(sm/h)
%     row_no= (i*(sm/h+1)-(sm/h)); %%where wants to insert
%     W(1:row_no-1,:) = W_prev(1:row_no-1,:) ;
%     tp =W_prev(row_no:end,:);
%     W(row_no,:)=b;
%     W(row_no+1:end+1,:) =tp;
%     W_prev = W;
%     
%     row_no=(i*(sm/h+1));
%     W(1:row_no-1,:) = W_prev(1:row_no-1,:) ;
%     tp =W_prev(row_no:end,:);
%     W(row_no,:)=b;
%     W(row_no+1:end+1,:) =tp;
%     W_prev = W;   
%     end
%     
%     W = W_prev;
%     
%     %W = vec2mat(W,(sm/h+1));
%     W = reshape(W,[sm/h+1,sm/h+1]);
%     W = W';
%     %W = flipud(W);
%     
% %     %Grid Matrix
%      x = [0:h:sm];
%      y = [0:h:sm];
%      [x,y] = meshgrid(x,y);
% 
%     %Plot
% %    if (time_intg_type ==1)
%        if (mod(ts,50)==0)
%     fig1 = figure (1); clf;               
%     hold on;        
%     temp = [num2str(t),' second - Implicit']; title(temp);     
%     
%     contourf(x,y,W,50,'LineColor','none')
%     %contourf(x,y,W)
%     
%     colorbar;
%     caxis([-1e-6 1e-6])
% %     
% %             filename =  ['graphs_u/response_at_' num2str(ts) '_step_' num2str(t) 's.pdf'];
% %             print(fig1, '-r90', '-dpdf', filename)
% %     
%     pause (1)
%        end 
%        
%        
%        
%     
% %     else
% %     if (mod(ts,5000)==0)
% %     fig1 = figure (1); clf;               
% %     hold on;        
% %     temp = [num2str(t),' second - Explicit']; title(temp);     
% %     
% %     contourf(x,y,W,50,'LineColor','none')
% %     %contourf(x,y,W)
% %     colorbar;
% %     %caxis([-15e-7 15e-7])
% %         
% % %             filename =  ['graphs_u_explicit/response_at_' num2str(ts) '_step_' num2str(t) 's.pdf'];
% % %             print(fig1, '-r90', '-dpdf', filename)
% %     
% %     pause (1)
% %     end
%     
%     end
end