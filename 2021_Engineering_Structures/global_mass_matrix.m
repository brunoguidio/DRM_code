function FGM = global_mass_matrix (n)
p = 2500;
A = 0.25*0.4;
I = 0.25*(0.4^3)/12;
global L;
h = L/n;


%Mass Matrix
%4x4
for i = 1:4
    for j = 1:4
        
        M3_LE_cmp = @(x)(sh_3(i,x).*sh_3(j,x))*(h/2)*p*A;
        
        M3_LE(i,j) = quad(M3_LE_cmp,-1,1);
    end
end

%3x3
for i = 1:3
    for j = 1:3
        
        M2_LE_cmp = @(x)(sh_2(i,x).*sh_2(j,x))*(h/2)*p*I;
        
        M2_LE(i,j) = quad(M2_LE_cmp,-1,1);
    end
end

%ELEMENT Mass Matrix
EMM = sparse(zeros (7,7));
EMM (1:4,1:4) = M3_LE;
EMM (5:7,5:7) = M2_LE;

%GLOBAL MASS MATRIX
for nd = 1:n
    GM = sparse(zeros (5*n+2,5*n+2));
    
    GM (3*nd-2:3*nd+1,3*nd-2:3*nd+1)= EMM (1:4,1:4);
    GM (3*n+2*nd:3*n+2*nd+2,3*n+2*nd:3*n+2*nd+2) = EMM (5:7,5:7);
    
    if nd == 1
        FGM = GM;
    else
        FGM = GM + GM_prev;
    end
    
    GM_prev = FGM;
    
end

%Boundary Conditions
FGM(:,3*n+1)=[]; FGM(3*n+1,:)=[]; FGM(:,2*n+1)=[]; FGM(2*n+1,:)=[];
FGM(:,n+1)=[]; FGM(n+1,:)=[]; FGM(:,1)=[]; FGM(1,:)=[];