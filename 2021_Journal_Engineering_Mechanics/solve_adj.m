function [l_bottom] = solve_adj(u)

global sm; global h; global rho; global G;
global T; global freq; global ne; global nn; global nb;
global time_intg_type; global dt; global beta; global gamma; 
global s_loc; global Ntimestep;
global GK; global GM; global Keff; global GM_inv; global um;

GF = sparse(zeros(nn-nb,1));
l_curt = zeros(length(GF),1);
ldot_curt = zeros(length(GF),1);
lddot_curt = GM_inv * GF;

for ts = Ntimestep:-1:1
    
    t = ts*dt; 
    l_prev = l_curt;
    ldot_prev = ldot_curt;
    lddot_prev = lddot_curt;
    
    %Global Force Vector
    for i = 1:length(s_loc)    
    GF(s_loc(i),1) = 2*(u(i,ts) - um(i,ts)) ;
    end

    %Build the RHS load vector
    if (time_intg_type ==1)
    %Implicit time integration
        RHS = GF - GK*(l_prev - ldot_prev*dt + 0.25*lddot_prev*dt^2);
        lddot_curt = Keff\RHS;
    else       
    %explicit time integration.
        RHS = GF - GK * (l_prev);
        lddot_curt = GM_inv * RHS;
    end
    
    %Updating u and udot at the current (i)-th time step
    l_curt = l_prev - ldot_prev * dt + ...
        0.5 * (lddot_prev * (1.0-2.0*beta) + lddot_curt* 2.0*beta)*dt^2;
    
    ldot_curt = ldot_prev - ...
        (lddot_prev * (1.0-gamma) + lddot_curt* gamma)*dt;
    
    l_curt_history(:,ts) = l_curt;
    
    %lambda all bottom nodes
    l_bottom(:,ts) = l_curt(1:(sm/h+1),1);
    
end