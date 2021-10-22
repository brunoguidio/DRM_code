%for i = GN 
clf; 
mean_error = mean(error_mat_each_it,2); % 2 means mean of a row, % 1 means mean of a column
%end

smallest_error = min(error_mat_each_it,[],2); % 2 means mean of a row, % 1 means mean of a column


fig = figure(1);
semilogy([1:GN+1],mean_error); grid on
    ylabel('mean value of $\mathcal{L}$','Interpreter','Latex',...
        'FontSize',12);
    xlabel('GA iterations (generations)','Interpreter','Latex',...
        'FontSize',12);
    
    set(gcf, 'Position', [200,200,600,400])
    filename =  ['mean_error_case2.png'];   
%     
    print(fig, '-r90', '-dpng', filename); 


fig = figure(2);
semilogy([1:GN+1],smallest_error); grid on
    ylabel('$\mathcal{L}$ for the best fit individual','Interpreter','Latex',...
        'FontSize',12);
    xlabel('GA iterations (generations)','Interpreter','Latex',...
        'FontSize',12);
    ylim ([10^-9 10^-5])
    
    set(gcf, 'Position', [200,200,600,400])
    filename =  ['smallest_error_case2.png'];   
%     
    print(fig, '-r90', '-dpng', filename); 





%%
[a1,a2,a3] = size(P_population_history);
for i = 1:a2
    for j = 1:a3
        P1(i,j) = P_population_history(1,i,j);
    end
end

fig_handle = figure(3);
bar3(P1);
 ylabel('GA iterations','Interpreter','Latex','FontSize',12);
 xlabel('Individuals','Interpreter','Latex','FontSize',12);
 zlabel('$P$ [$\textrm{N}/\textrm{m}$]','Interpreter','Latex','FontSize',12); 
 view([-75 09])
 set(gcf, 'Position', [200,200,600,400])
       filename =  ['P_history.png'];   
       print(fig_handle, '-r190', '-dpng', filename)

%%

[a1,a2,a3] = size(freq_population_history);
for i = 1:a2
    for j = 1:a3
        freq1(i,j) = freq_population_history(1,i,j);
    end
end

fig_handle = figure(4);
bar3(freq1);
 ylabel('GA iterations','Interpreter','Latex','FontSize',12);
 xlabel('Individuals','Interpreter','Latex','FontSize',12);
 zlabel('$f$ [Hz]','Interpreter','Latex','FontSize',12); 
 view([-75 09])
 %set(gcf, 'Position', [200,200,600,400])
       filename =  ['freq_history.png'];   
       print(fig_handle, '-r190', '-dpng', filename)


[a1,a2,a3] = size(E1_population_history);
for i = 1:a2
    for j = 1:a3
        E1(i,j) = E1_population_history(1,i,j)*10^10;
    end
end

fig_handle = figure(5);
bar3(E1);
 ylabel('GA iterations','Interpreter','Latex','FontSize',12);
 xlabel('Individuals','Interpreter','Latex','FontSize',12);
 zlabel('$E_1$ [Pa]','Interpreter','Latex','FontSize',12); 
 view([-75 09])
 %set(gcf, 'Position', [200,200,600,400])
       filename =  ['E1_history.png'];   
       print(fig_handle, '-r190', '-dpng', filename)

[a1,a2,a3] = size(E2_population_history);
for i = 1:a2
    for j = 1:a3
        E2(i,j) = E2_population_history(1,i,j)*10^10;
    end
end

fig_handle = figure(6);
bar3(E2);
 ylabel('GA iterations','Interpreter','Latex','FontSize',12);
 xlabel('Individuals','Interpreter','Latex','FontSize',12);
 zlabel('$E_2$ [Pa]','Interpreter','Latex','FontSize',12); 
 view([-75 09])
 %set(gcf, 'Position', [200,200,600,400])
       filename =  ['E2_history.png'];   
       print(fig_handle, '-r190', '-dpng', filename)
       
[a1,a2,a3] = size(E3_population_history);
for i = 1:a2
    for j = 1:a3
        E3(i,j) = E3_population_history(1,i,j)*10^10;
    end
end

fig_handle = figure(7);
bar3(E3);
 ylabel('GA iterations','Interpreter','Latex','FontSize',12);
 xlabel('Individuals','Interpreter','Latex','FontSize',12);
 zlabel('$E_3$ [Pa]','Interpreter','Latex','FontSize',12); 
 view([-75 09])
% set(gcf, 'Position', [200,200,600,400])
       filename =  ['E3_history.png'];   
       print(fig_handle, '-r190', '-dpng', filename)

% [a1,a2,a3] = size(E4_population_history);
% for i = 1:a2
%     for j = 1:a3
%         E4(i,j) = E4_population_history(1,i,j);
%     end
% end
% 
% figure(8);
% bar3(E4);
% 
% [a1,a2,a3] = size(E5_population_history);
% for i = 1:a2
%     for j = 1:a3
%         E5(i,j) = E5_population_history(1,i,j);
%     end
% end
% 
% figure(9);
% bar3(E5);
% 
% [a1,a2,a3] = size(E6_population_history);
% for i = 1:a2
%     for j = 1:a3
%         E6(i,j) = E6_population_history(1,i,j);
%     end
% end
% 
% figure(10);
% bar3(E6);
% 
% [a1,a2,a3] = size(E7_population_history);
% for i = 1:a2
%     for j = 1:a3
%         E7(i,j) = E7_population_history(1,i,j);
%     end
% end
% 
% figure(11);
% bar3(E7);
% 
% [a1,a2,a3] = size(E8_population_history);
% for i = 1:a2
%     for j = 1:a3
%         E8(i,j) = E8_population_history(1,i,j);
%     end
% end
% 
% figure(12);
% bar3(E8);
% 
% [a1,a2,a3] = size(E9_population_history);
% for i = 1:a2
%     for j = 1:a3
%         E9(i,j) = E9_population_history(1,i,j);
%     end
% end
% 
% figure(13);
% bar3(E9);

% [a1,a2,a3] = size(E10_population_history);
% for i = 1:a2
%     for j = 1:a3
%         E10(i,j) = E10_population_history(1,i,j);
%     end
% end
% 
% figure(14);
% bar3(E10);



