function [Rr] = row_regularization(p,R,a,b1,b2,b3,b4,c,nM,nRM)
%global a...
 

% a = [1,nM,nM*nRM+1,nM*(nRM+1)];
% b1 = [2:(nM-1)]; b2 = [nM+1:nM:nM*nRM+1-nM]; b3 = [nM*2:nM:nM*(nRM+1)-nM]; b4 = [nM*nRM+2:nM*(nRM+1)-1];
% c = [7,8,9,12,13,14,17,18,19];



Rr(1,1:nM*(nRM+1)) = 0;
   Lia = ismember(p,a);
   if Lia == 1;
        if p == 1; %row1
        Rr(1,[p,p+1,p+nM,p+nM+1]) = [R(1) R(2) R(4) R(3)]; 

        elseif p == nM;%row5
        Rr(1,[p-1,p,p+nM-1,p+nM]) = [R(2) R(1) R(3) R(4)];

        elseif p == nM*nRM+1; %row21
        Rr(1,[p-nM,p-nM+1,p,p+1]) = [R(4) R(3) R(1) R(2)];

        elseif p == nM*(nRM+1); %row25
        Rr(1,[p-nM-1,p-nM,p-1,p]) = [R(3) R(4) R(2) R(1)];
        end
        
   else
       Lia = ismember(p,b1);
       if Lia == 1;
       Rr(1,[p-1, p, p+1, p+nM-1, p+nM, p+nM+1]) = [R(2) 2*R(1) R(2) R(3) 2*R(4) R(3)];
       
       else
           Lia = ismember(p,b2);
           if Lia == 1;
           Rr(1,[p-nM, p-nM+1, p, p+1,p+nM, p+nM+1]) = [R(4) R(3) 2*R(1) 2*R(2) R(4) R(3)];
           
           else 
               Lia = ismember(p,b3);
               if Lia == 1;
                   Rr(1,[p-nM-1, p-nM, p-1, p, p+nM-1,p+nM]) = [R(3) R(4) 2*R(2) 2*R(1) R(3) R(4)];
               
               else
                   Lia =ismember(p,b4);
                   if Lia == 1;
                       Rr(1,[p-nM-1, p-nM, p-nM+1,p-1, p,p+1]) = [R(3) 2*R(4) R(3) R(2) 2*R(1) R(2)];
                   
                   else
                       Lia = ismember(p,c);
                       if Lia == 1;
                       Rr(1,[p-nM-1, p-nM, p-nM+1, p-1, p, p+1, p+nM-1, p+nM, p+nM+1]) = [R(3) 2*R(4) R(3) 2*R(2) 4*R(1) 2*R(2) R(3) 2*R(4) R(3)];
                       else
                           Rr = 0;
                       end
                   end
               end
           end
       end
   end    
