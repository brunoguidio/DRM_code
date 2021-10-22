function [um_noised] = Function_add_noise_to_um (um,noise_level)

[row col] = size(um);

for i = 1:row
       um_noised(i,:) = um(i,:) + randn(size(um(i,:)))*max(um(i,:))*noise_level/100;
end