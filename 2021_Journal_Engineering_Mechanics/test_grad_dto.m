
hold on
plot(1:1500,um(15,:))
plot(1:1500,u(15,:),'--r')

hold on
plot(1:1500,l_bottom(31,:))
plot(1:1500,l_RBL(31,:),'r--')

grad_otd = -h*dt*l_bottom;
grad_dto = h*l_RBL;

hold on
plot(1:1500,grad_otd(31,:))
plot(1:1500,grad_dto(31,:),'r--')
%%
figure(1)
l_bottom_n = grad_otd./abs(max(grad_otd)); 
l_RBL_n = grad_dto./abs(max(grad_dto));

hold on
plot([1:1500],l_bottom_n(31,:))
plot([1:1500],l_RBL_n(31,:),'--r')
legend('adj','adj_{dto}')

%%
figure(2)
l_bottom_2 = grad_otd./sqrt(dot(grad_otd,grad_otd));
l_RBL_n_2 = grad_dto./sqrt(dot(grad_dto,grad_dto));

hold on %Ntimestep = 1500;
plot([1:Ntimestep],l_bottom_2(31,:))
plot([1:Ntimestep],l_RBL_n_2(31,:),'--r')
legend('adj','adj_{dto}')