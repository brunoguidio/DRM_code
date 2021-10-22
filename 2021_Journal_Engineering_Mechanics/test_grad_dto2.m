hold on
plot(1:1500,grad(31,:))
plot(1:1500,grad_dto(31,:),'r--')

%%
grad_n = grad/abs(max(grad(:)));
grad_dto_n = grad_dto/abs(max(grad_dto(:)));

hold on
plot(1:1500,grad_n(31,:))
plot(1:1500,grad_dto_n(31,:),'r--')

%%
grad_n_2 = grad./sqrt(dot(grad,grad));
grad_dto_n_2 = grad_dto./sqrt(dot(grad_dto,grad_dto));


hold on
plot(1:1500,grad_n_2(31,:))
plot(1:1500,grad_dto_n_2(31,:),'r--')
