%Prueba para pscglobal
clc;clear all;close all;
fx =  @(x) -x(1)*x(2)*x(3);
hx = @(x) x(1)*x(2) + x(1)*x(3) + x(2)*x(3) -27;
x0 = [3;3;3];

%% Resultados:
[x, lambda, k] = pcsglobal(fx, hx, x0);

X = x(1)
Y = x(2)
Z = x(3)
