function y = fine_discretize(truecp,Ldt,Sdt)

%Ldt : larger dt
%Sdt : smaller dt

[rows cols] = size(truecp);

if (rem(Ldt,Sdt)==0)
    m = 0; 
    temp = Ldt/Sdt;
    for i = 1:rows-1
        for j = 1:temp
            m = m+1;
            y(m,1) =  truecp(i) + (truecp(i+1)-truecp(i))/temp*(j-1); 
        end 
    end
else
    fprintf('larger dt is not a multiple of smaller dt::::');
    temp = input('Enter Ctrl + C;')
end 

       
return