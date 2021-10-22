function er = error_movingload(y)
    
global ne;
global ns;
global u_true;
global iteration_so_far; 
global error_mat_each_it;
global GN; global PS; global L; 
global P_population_history; global freq_population_history; 
global E1_population_history;global E2_population_history;
global E3_population_history;%global E4_population_history;
% global E5_population_history;global E6_population_history;
% global E7_population_history;global E8_population_history;
% global E9_population_history;%global E10_population_history;
global N_GS_Load; 

Pg = [y(1:N_GS_Load)];
freqg = [y(N_GS_Load+1:N_GS_Load*2)];
E1g = [y(N_GS_Load*2+1)];
E2g = [y(N_GS_Load*2+2)];
E3g = [y(N_GS_Load*2+3)];
% E4g = [y(N_GS_Load*2+4)];
% E5g = [y(N_GS_Load*2+5)];
% E6g = [y(N_GS_Load*2+6)];
% E7g = [y(N_GS_Load*2+7)];
% E8g = [y(N_GS_Load*2+9)];
% E9g = [y(N_GS_Load*2+9)];
%E10g = [y(N_GS_Load*2+10)];

u_guess = CB_ML_5 (Pg,freqg,E1g,E2g,E3g,ne,ns);  
       
    %CJ
    [row col] = size(u_true); 
    z = zeros(row,col);
    for k = 1:row
        for j = 1:col
            z(k,j) = (u_true(k,j)-u_guess(k,j))^2;
        end 
    end 
    
    er = sum(sum(z))
 
    iteration_so_far =  iteration_so_far + 1
     
    row = ceil(iteration_so_far/PS);
    col = mod(iteration_so_far,PS);%its remainder  
    if (col==0) 
        col = PS; 
    end  

    for i = 1:N_GS_Load
        P_population_history(i,row,col) = Pg(i);
        freq_population_history(i,row,col) = freqg(i);
        E1_population_history(i,row,col) = E1g;
        E2_population_history(i,row,col) = E2g;
        E3_population_history(i,row,col) = E3g;
%         E4_population_history(i,row,col) = E4g;
%         E5_population_history(i,row,col) = E5g;
%         E6_population_history(i,row,col) = E6g;
%         E7_population_history(i,row,col) = E7g;
%         E8_population_history(i,row,col) = E8g;
%         E9_population_history(i,row,col) = E9g;
%         E10_population_history(i,row,col) = E10g;
    end 

    error_mat_each_it(row,col) = er;
    
    %fprintf('error_fit \n');
end 