function [l_RBL] = solve_adj_dto(u)

global h; 
global nn;
global dt; 
global s_loc; global Ntimestep; 
global GK; global GM; global um;
%global ln; global bn; global rn; global wm; global hm; global Invert_only;
global sm; %new 1/17/20


a0 = 4/dt^2; a1 = 2/dt; a2 = 4/dt; a3 = 1; a4 = 1; a5=0;  
GC = sparse(zeros(size(GM))); %new 1/17/20

Keff_dto = sparse(a0*GM+a1*GC+GK);
Keff_t = Keff_dto';
GM_t = GM';
GC_t = GC';

%Bd = sparse(zeros(nn,1));
[row col] = size(GM); %new
Bd = sparse(zeros(row,1)); %new

for ts = Ntimestep:-1:1
    
 %   t = ts*dt; 

    if ts == Ntimestep
%     lddot_curt = (zeros(nn,1));
%     ldot_curt = (zeros(nn,1));
    lddot_curt = (zeros(row,1)); %new
    ldot_curt = (zeros(row,1)); %new
    
    Bd(s_loc,1) = -2*(u(:,ts) - um(:,ts)) ;
    
    upL = dt*Bd + a1*ldot_curt +a0*lddot_curt;
  
    l_curt = Keff_t\upL;
    
    
%     elseif ts == 1
%         
%     upI = (a3*GM_t + a5*GC_t)*l_curt - a5*ldot_curt - a3*lddot_curt;
%     lddot_prev = GM_t\upI;
%     
%     ldot_prev= -GC_t*lddot_prev+(a2*GM_t + a4*GC_t)*l_curt - a4*ldot_curt - a2*lddot_curt;
%         
%     Bd(s_loc,1) = (u(:,ts) - um(:,ts)) ;
%     l_prev= -GK'*lddot_prev+(a0*GM_t + a1*GC_t)*l_curt - a1*ldot_curt - a0*lddot_curt+dt*Bd;  
%         
%     l_curt = l_prev;
%     ldot_curt = ldot_prev;
%     lddot_curt = lddot_prev;
    
    else
    lddot_prev = (a3*GM_t + a5*GC_t)*l_curt - a5*ldot_curt - a3*lddot_curt;
    
    ldot_prev = (a2*GM_t + a4*GC_t)*l_curt - a4*ldot_curt - a2*lddot_curt;
    
    Bd(s_loc,1) = -2*(u(:,ts) - um(:,ts));
    upL = dt*Bd + a1*ldot_prev + a0*lddot_prev + ...
          (a0*GM_t + a1*GC_t)*l_curt - a1*ldot_curt - a0*lddot_curt;
    
    l_prev = Keff_t\upL;
    
    
    l_curt = l_prev;
    ldot_curt = ldot_prev;
    lddot_curt = lddot_prev;
    end
    
   %lambda all bottom nodes
    l_RBL(:,ts) = l_curt(1:(sm/h+1),1);
    
    %l_RBL(:,ts) = l_curt([ln,bn,rn],1);
% if (strcmp(Invert_only,'on,off,off'))
%     l_RBL(:,ts) = l_curt([ln,1],1);
%     
% elseif (strcmp(Invert_only,'off,on,off'))
%     l_RBL(:,ts) = l_curt(bn,1);
%     
% elseif (strcmp(Invert_only,'off,off,on'))
%     l_RBL(:,ts) = l_curt([(wm/h+1),rn],1);
% 
% elseif (strcmp(Invert_only,'on,on,off'))
%     l_RBL(:,ts) = l_curt([ln,bn],1);
%     
% elseif (strcmp(Invert_only,'off,on,on'))
%     l_RBL(:,ts) = l_curt([bn,rn],1);
% 
% elseif (strcmp(Invert_only,'on,off,on'))
%     l_RBL(:,ts) = l_curt([ln,1,(wm/h+1),rn],1);
%       
% elseif (strcmp(Invert_only,'on,on,on'))
%     l_RBL(:,ts) = l_curt([ln,bn,rn],1);
% end
    
end