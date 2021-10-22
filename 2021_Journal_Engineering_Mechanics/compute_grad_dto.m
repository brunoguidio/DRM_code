function [grad_dto] = compute_grad_dto(F_g,l_RBL)

global h; global dt; global nn; 
global dt; global R_factor;global R_factor_intensity; 
global Fg_vector; global Reg_vector;
dx = h;

fun11 = @(n,e) dt/dx * ((n-1)/4).*((n-1)/4) + dx/dt * ((e-1)/4).*((e-1)/4);
value11 = dblquad(fun11,-1,1,-1,1);

fun12 = @(n,e) dt/dx * ((n-1)/4).*((1-n)/4) + dx/dt * ((e-1)/4).*(1/4*(-e-1));
value12 = dblquad(fun12,-1,1,-1,1);

fun13 = @(n,e) dt/dx * ((n-1)/4).*((n+1)/4) + dx/dt * ((e-1)/4).*((e+1)/4);
value13 = dblquad(fun13,-1,1,-1,1);

fun14 = @(n,e) dt/dx * ((n-1)/4).*(1/4*(-n-1)) + dx/dt * ((e-1)/4).*((1-e)/4);
value14 = dblquad(fun14,-1,1,-1,1);

% fun21 = @(n,e) dt/dx * ((1-n)/4).*((n-1)/4) + dx/dt * (1/4*(-e-1)).*((e-1)/4);
% value21 = dblquad(fun21,-1,1,-1,1);
% 
% fun22 = @(n,e) dt/dx * ((1-n)/4).*((1-n)/4) + dx/dt * (1/4*(-e-1)).*(1/4*(-e-1));
% value22 = dblquad(fun22,-1,1,-1,1);
% 
% fun23 = @(n,e) dt/dx * ((1-n)/4).*((n+1)/4) + dx/dt * (1/4*(-e-1)).*((e+1)/4);
% value23 = dblquad(fun23,-1,1,-1,1);
% 
% fun24 = @(n,e) dt/dx * ((1-n)/4).*(1/4*(-n-1)) + dx/dt * (1/4*(-e-1)).*((1-e)/4);
% value24 = dblquad(fun24,-1,1,-1,1);
% 
% fun31 = @(n,e) dt/dx * ((n+1)/4).*((n-1)/4) + dx/dt * ((e+1)/4).*((e-1)/4);
% value31 = dblquad(fun31,-1,1,-1,1);
% 
% fun32 = @(n,e) dt/dx * ((n+1)/4).*((1-n)/4) + dx/dt * ((e+1)/4).*(1/4*(-e-1));
% value32 = dblquad(fun32,-1,1,-1,1);
% 
% fun33 = @(n,e) dt/dx * ((n+1)/4).*((n+1)/4) + dx/dt * ((e+1)/4).*((e+1)/4);
% value33 = dblquad(fun33,-1,1,-1,1);
% 
% fun34 = @(n,e) dt/dx * ((n+1)/4).*(1/4*(-n-1)) + dx/dt * ((e+1)/4).*((1-e)/4);
% value34 = dblquad(fun34,-1,1,-1,1);
% 
% fun41 = @(n,e) dt/dx * (1/4*(-n-1)).*((n-1)/4) + dx/dt * ((1-e)/4).*((e-1)/4);
% value41 = dblquad(fun41,-1,1,-1,1);
% 
% fun42 = @(n,e) dt/dx * (1/4*(-n-1)).*((1-n)/4) + dx/dt * ((1-e)/4).*(1/4*(-e-1));
% value42 = dblquad(fun42,-1,1,-1,1);
% 
% fun43 = @(n,e) dt/dx * (1/4*(-n-1)).*((n+1)/4) + dx/dt * ((1-e)/4).*((e+1)/4);
% value43 = dblquad(fun43,-1,1,-1,1);
% 
% fun44 = @(n,e) dt/dx * (1/4*(-n-1)).*(1/4*(-n-1)) + dx/dt * ((1-e)/4).*((1-e)/4);
% value44 = dblquad(fun44,-1,1,-1,1);

% R = [value11 value12 value13 value14;
%      value12 value22 value23 value24;
%      value13 value23 value33 value34;
%      value14 value24 value34 value44];

R = [value11 value12 value13 value14];

%Mapping Matrix
[row col] = size(F_g);
nM = col;  %number of mapping matrix
nRM = row-1;   %number of rows for each mapping matrix

%%
a = [1,nM,nM*nRM+1,nM*(nRM+1)];
b1 = [2:(nM-1)]; b2 = [nM+1:nM:nM*nRM+1-nM]; 
b3 = [nM*2:nM:nM*(nRM+1)-nM]; b4 = [nM*nRM+2:nM*(nRM+1)-1];
c = [1:nM*(nRM+1)];

for i = length(c):-1:1;
Lia = ismember(c(i),a);
if Lia == 1
    c(i) = [];
end
end

for i = length(c):-1:1;
Lia = ismember(c(i),b1);
if Lia == 1
    c(i) = [];
end
end

for i = length(c):-1:1;
Lia = ismember(c(i),b2);
if Lia == 1
    c(i) = [];
end
end

for i = length(c):-1:1;
Lia = ismember(c(i),b3);
if Lia == 1
    c(i) = [];
end
end

for i = length(c):-1:1;
Lia = ismember(c(i),b4);
if Lia == 1
    c(i) = [];
end
end

%%
Fg_vector = reshape(F_g',nM*(nRM+1),1);
for p = 1:nM*(nRM+1)
    
    Rr = row_regularization(p,R,a,b1,b2,b3,b4,c,nM,nRM);
    Reg_vector(p,1) = Rr*Fg_vector;
    
    %GR2(p,:) = Rr;
    
end

Reg = reshape(Reg_vector,nM,nRM+1)';

%%
R_factor = R_factor_intensity * sqrt(dot(l_RBL',l_RBL'))/...
    sqrt(dot(Reg,Reg)+10^(-50));

%%
grad_dto = h*l_RBL + R_factor* Reg';