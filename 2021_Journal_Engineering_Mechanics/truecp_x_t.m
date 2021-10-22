function [truecp] = truecp_x_t (truecp)
    %i = 2:21; j = 22:41; k = 42:61;
global sm; global h; global Ntimestep; 

    dk = 30;
for i = 2:sm/h+1
    if i <= (sm/h)/3
        dk(i) = dk(i-1)-1;
        
    elseif i <= (sm/h)*(2/3)
        dk(i) = dk(i-1)+1;
        
    else
        dk(i) = dk(i-1)-1;
    end
    
end
dk(1) = 0;

    for i = 2:sm/h+1
    truecp(sum(dk(1:i))+1:end,i) = truecp(1:end-(sum(dk(1:i))),1);
    end
    
%    x = [1:sm/h+1];
%    y = [1:Ntimestep];
%    [x,y] = meshgrid(x,y); 
%    contourf(x,y,truecp,50,'LineColor','none')
%    colorbar;
%    caxis([-1 2])