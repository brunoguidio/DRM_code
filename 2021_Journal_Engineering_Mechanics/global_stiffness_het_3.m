function [GK] = global_stiffness_het_3 (MM)

global sm; global h; global rho; global G;
global T; global freq; global ne; global nn; global nb;
global time_intg_type; global dt; global beta; global gamma; 
global s_loc; global Ntimestep;
global GM; global Keff; global GM_inv;


%Stiffness Matrix
Ke = [0.6667   -0.1667   -0.3333   -0.1667;
     -0.1667    0.6667   -0.1667   -0.3333;
     -0.3333   -0.1667    0.6667   -0.1667;
     -0.1667   -0.3333   -0.1667    0.6667];
Ke = Ke; %(h^2);
 
%Build Global Stiffness Matrix
GK = sparse(nn,nn);

%GK = zeros(nn,nn);
%Heterogeneous Domain
%Location 
L1 = [0,60;0,30];
L2 = [25,35;25,35];
L3 = [0,60;30,60];
%L4 = [10,20;20,30];
%L5 = [0,60;40,60];
%L6 = [30,40;40,50];

%% Part 1 and 2
Le1 = [L1(1,1)/h+1+(sm/h)*L1(2,1)/h:(L1(1,2)-h)/h+1+(sm/h)*L1(2,1)/h];
for i = 2:(L1(2,2)-L1(2,1))/h
  Le1(i,:) = Le1(i-1,:) + (sm/h);  
end
Le1 = Le1(:);

Le2 = [L2(1,1)/h+1+(sm/h)*L2(2,1)/h:(L2(1,2)-h)/h+1+(sm/h)*L2(2,1)/h];
for i = 2:(L2(2,2)-L2(2,1))/h
  Le2(i,:) = Le2(i-1,:) + (sm/h);  
end

for i = length(Le1):-1:1
    Lia = ismember(Le1(i),Le2);
    if Lia == 1
        Le1(i) = [];
    end
end
Le2 = Le2(:);

for k = 1:length(Le1)
    e = Le1(k);
    for i = 1:4
    for j = 1:4
        mi = MM(e,i); 
        mj = MM(e,j);
        GK(mi,mj) = Ke(i,j)*G(1)+ GK(mi,mj);
    end
    end
end 

for k = 1:length(Le2)
    e = Le2(k);
    for i = 1:4
    for j = 1:4
        mi = MM(e,i); 
        mj = MM(e,j);
        GK(mi,mj) = Ke(i,j)*G(2)+ GK(mi,mj);
    end
    end
end 
 %% Part 3 and 4
 Le3 = [L3(1,1)/h+1+(sm/h)*L3(2,1)/h:(L3(1,2)-h)/h+1+(sm/h)*L3(2,1)/h];
for i = 2:(L3(2,2)-L3(2,1))/h
  Le3(i,:) = Le3(i-1,:) + (sm/h);  
end
Le3 = Le3(:);

Le2 = [L2(1,1)/h+1+(sm/h)*L2(2,1)/h:(L2(1,2)-h)/h+1+(sm/h)*L2(2,1)/h];
for i = 2:(L2(2,2)-L2(2,1))/h
  Le2(i,:) = Le2(i-1,:) + (sm/h);  
end

for i = length(Le3):-1:1
    Lia = ismember(Le3(i),Le2);
    if Lia == 1
        Le3(i) = [];
    end
end
Le2 = Le2(:);

for k = 1:length(Le3)
    e = Le3(k);
    for i = 1:4
    for j = 1:4
        mi = MM(e,i); 
        mj = MM(e,j);
        GK(mi,mj) = Ke(i,j)*G(3)+ GK(mi,mj);
    end
    end
end 

% for k = 1:length(Le4)
%     e = Le4(k);
%     for i = 1:4
%     for j = 1:4
%         mi = MM(e,i); 
%         mj = MM(e,j);
%         GK(mi,mj) = Ke(i,j)*G(4)+ GK(mi,mj);
%     end
%     end
% end 

% % Part 5 and 6
% Le5 = [L5(1,1)/h+1+(sm/h)*L5(2,1)/h:(L5(1,2)-h)/h+1+(sm/h)*L5(2,1)/h];
% for i = 2:(L5(2,2)-L5(2,1))/h
%   Le5(i,:) = Le5(i-1,:) + (sm/h);  
% end
% Le5 = Le5(:);
% 
% Le6 = [L6(1,1)/h+1+(sm/h)*L6(2,1)/h:(L6(1,2)-h)/h+1+(sm/h)*L6(2,1)/h];
% for i = 2:(L6(2,2)-L6(2,1))/h
%   Le6(i,:) = Le6(i-1,:) + (sm/h);  
% end
% 
% for i = length(Le5):-1:1
%     Lia = ismember(Le5(i),Le6);
%     if Lia == 1
%         Le5(i) = [];
%     end
% end
% Le6 = Le6(:);
% 
% for k = 1:length(Le5)
%     e = Le5(k);
%     for i = 1:4
%     for j = 1:4
%         mi = MM(e,i); 
%         mj = MM(e,j);
%         GK(mi,mj) = Ke(i,j)*G(5)+ GK(mi,mj);
%     end
%     end
% end 
% 
% for k = 1:length(Le6)
%     e = Le6(k);
%     for i = 1:4
%     for j = 1:4
%         mi = MM(e,i); 
%         mj = MM(e,j);
%         GK(mi,mj) = Ke(i,j)*G(6)+ GK(mi,mj);
%     end
%     end
% end 
%  