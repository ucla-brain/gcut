function [fitness] = Gassi_com_fit_beta(angle, order)
%This function is used to compute fitness of a neurite for a given soma.

mean = 1.13;
a = 1;

std = 0.69;
fun =@(x) 2*a*exp((-(x-mean).^2)/(2*std.^2));
fitness = integral(fun, angle, pi);
fitness = fitness * (1+log(1+4/order));


end


