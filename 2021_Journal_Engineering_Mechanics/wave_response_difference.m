   W2 = W_target(1:2:121,1:2:121);

   error = abs(W2-W_guess).^2;
   error_norm_history_north = sum(error)/sum(abs(W2).^2)*100;
   
   
   
     figure_F = figure (2);
 %   subplot(1,2,2)
    hold on
%     temp = ['Case 2']; 
%     title(temp); 
    contourf(x,y,W2-W_guess,50,'LineColor','none') 
    colorbar;
    caxis([-4e-6 4e-6])
    
    ylabel('$y$ [m]','Interpreter','Latex',...
        'FontSize',12);
  
    xlabel('$x$ [m]','Interpreter','Latex',...
        'FontSize',12);
    
   filename =  ['Difference_wave_responses_case3.pdf'];   
%     
   print(figure_F, '-r100', '-dpdf', filename); 