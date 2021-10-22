function [GK, GM, GF_i] = FEM_matrices ( )

global sm; global h; global rho; global G;
global T; global freq; global ne; global nn; global nb;
global time_intg_type; global dt; global beta; global gamma; 
global s_loc; global Ntimestep;
global Keff; global GM_inv; global MM; global BC; global Soil_profile;

%Mapping Matrix
nM = (sm/h)+1;  %number of mapping matrix
nRM = (sm/h);   %number of rows for each mapping matrix

MM = sparse(zeros((nM-1)*nRM,4)); %Sparse
MM(1,1:4) = [1, 2, nM+2, nM+1]; %[nM+1, nM+2, 2, 1];
for i = 2:nRM
    MM(i,1:4) = MM(i-1,1:4)+1;
end

for i = 2:nM-1
    MM(nRM*(i-1)+1:nRM*i,1:4) = MM(nRM*(i-2)+1:nRM*(i-1),1:4)+(sm/h+1);
end

GM = sparse(nn,nn);
GM_inv = sparse(nn,nn);
Keff = sparse(nn,nn);

%%%%%%%%%%  Global K  %%%%%%%%%%
if (strcmp(Soil_profile,'Homogeneous'))
    GK = global_stiffness_homo(MM); 
elseif (strcmp(Soil_profile,'Heterogeneous_3'))
    GK = global_stiffness_het_3(MM);
elseif (strcmp(Soil_profile,'Heterogeneous_6'))
    GK = global_stiffness_het_6(MM);
end
%GK = global_stiffness_homo(MM);

    
%Mass Matrix
Me =  [0.4444    0.2222    0.1111    0.2222;
       0.2222    0.4444    0.2222    0.1111;
       0.1111    0.2222    0.4444    0.2222;
       0.2222    0.1111    0.2222    0.4444];
[row, col] = size(Me);   

if (time_intg_type ==1) 
    Me = rho * ((h^2)/4) * Me; 
    
else %For Explicit    
Me_lumping = sparse(zeros(row,col));
for i = 1:row
    for j = 1:col
            Me_lumping(i,i) = Me_lumping(i,i) + Me(i,j);
    end
end
Me = rho * ((h^2)/4) * Me_lumping;   
end
  
%Build Global Mass Matrix
for e = 1:ne
    for i = 1:4
    for j = 1:4
        mi = MM(e,i); 
        mj = MM(e,j); 
        GM(mi,mj) = Me(i,j)+ GM(mi,mj);
    end
    end
end 

%Global Mass Matrix - Inverse
[row, col] = size(GM);

if (time_intg_type ==1)
    %GM_inv = inv(GM);
    GM_inv = zeros(row,col);
else
    for i = 1:row
    GM_inv(i,i)= 1/GM(i,i);
    end
end

%%%% GF_i %%%
%Global Force Vector
%t = 0;
GF_i = sparse(zeros(nn,1));

F = [1; 1; 0; 0]* h/2; 
for e = 1:(sm/h)
    for i = 1:4
    mi = MM(e,i); 
    GF_i(mi,1) = F(i,1)+ GF_i(mi,1);
    end
end


%%%%%%%%%%%%%% Boundary Conditions %%%%%%%%%%%%%%
BC = [1:nn]';
for i = (sm/h):-1:2
    BC(i*(sm/h+1))=[];
    BC(i*(sm/h+1)-(sm/h))=[];
end

GK = GK(BC,BC);
GM = GM(BC,BC);
GM_inv = GM_inv(BC,BC);
GF_i = GF_i(BC,1);
