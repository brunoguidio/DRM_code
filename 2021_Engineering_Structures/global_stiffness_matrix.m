function FGK = global_stiffness_matrix (n,E_1,E_2,E_3)
A = 0.25*0.4;
I = 0.25*(0.4^3)/12;
%G = E/2.5;
Ks = 5/6;
global L;
h = L/n;

E = [E_1;E_2;E_3]*10^10;
he = length(E);

%Stiffness Matrix
%3x3
for i = 1:3
    for j = 1:3
        
        K2_LE_cmp = @(x)(sh_2(i,x).*sh_2(j,x))*h/2;
        temp1 = A*Ks*quad(K2_LE_cmp,-1,1); %%
        
        K2_LE_cmp = @(x)(sh_2_1(i,x).*sh_2_1(j,x))*2/h;
        temp2 = 2.5*I* quad(K2_LE_cmp,-1,1); %%
        
        K2_LE(i,j) = temp1 + temp2;
        
    end
end

%4x4
for i = 1:4
    for j = 1:4
        
        K3_LE_cmp = @(x)(sh_3_1(i,x).*sh_3_1(j,x))*(2/h)*A*Ks;%%
        
        K3_LE(i,j) = quad(K3_LE_cmp,-1,1,1);
    end
end

%3x4
for i = 1:3
    for j = 1:4
        
        K23_LE_cmp = @(x)(sh_2(i,x).*sh_3_1(j,x))*A*Ks; %%
        
        K23_LE(i,j) = -quad(K23_LE_cmp,-1,1);
    end
end

%4x3
for i = 1:4
    for j = 1:3
        
        K32_LE_cmp = @(x)(sh_3_1(i,x).*sh_2(j,x))*A*Ks;
        
        K32_LE(i,j) = -quad(K32_LE_cmp,-1,1);
    end
end

% %ELEMENT Stifness Matrix
% ESM = sparse(zeros (7,7));
% ESM (1:4,1:4) = K3_LE;
% ESM (1:4,5:7) = K32_LE;
% ESM (5:7,1:4) = K23_LE;
% ESM (5:7,5:7) = K2_LE;

%GLOBAL STIFNEES MATRIX
nn = n/he;
for i = 1:he
for nd = nn*(i-1)+1:nn*i
    GK = sparse(zeros (5*n+2,5*n+2));
    %GK = zeros (5*n+2,5*n+2);   
    
    GK (3*nd-2:3*nd+1,3*nd-2:3*nd+1)= K3_LE*E(i)/2.5;
    GK (3*nd-2:3*nd+1,3*n+2*nd:3*n+2*nd+2) = K32_LE*E(i)/2.5;
    GK (3*n+2*nd:3*n+2*nd+2,3*nd-2:3*nd+1) = K23_LE*E(i)/2.5;
    GK (3*n+2*nd:3*n+2*nd+2,3*n+2*nd:3*n+2*nd+2) = K2_LE*E(i)/2.5;  
    
    if nd == 1
        FGK = GK;
    else
        FGK = GK + GK_prev;
    end
    GK_prev = FGK;
    
end
end

%Boundary Conditions
FGK(:,3*n+1)=[]; FGK(3*n+1,:)=[]; FGK(:,2*n+1)=[]; FGK(2*n+1,:)=[];
FGK(:,n+1)=[]; FGK(n+1,:)=[]; FGK(:,1)=[]; FGK(1,:)=[];