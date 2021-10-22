
function y = sh_2(i,x)
%y = (x-1/2).*(x+1/2)*2;

switch(i)
    
    case(1) 
        y = x.*(x-1)/2;
    case(2) 
        y = (1-x.*x);
    case(3) 
        y = x.*(x+1)/2;

end 

return