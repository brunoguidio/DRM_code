function [L] = Function_errorf_dto(um,u)

global dt; global R_factor; global Fg_vector; global Reg_vector;

[row,col] = size(um);
um_vector = reshape(um,row*col,1);
u_vector = reshape(u,row*col,1);

L1 = ((u_vector-um_vector)'*dt)*(u_vector-um_vector);

%Fg_vector = reshape(F_g',nM*(nRM+1),1);


L2 = Fg_vector'*Reg_vector;

L = L1+ (R_factor/2) *L2;