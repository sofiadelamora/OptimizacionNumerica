function [g,h] = hesfera(x)
% Funci�n de restricciones del problema de np puntos en la esfera unitaria
% de dimensi�n tres.
%
% Optimizaci�n Num�rica
% ITAM
% 20 de octubre de 2020

g = []; % uso de fmincon.m en Matlab

n = length(x);
np = floor(n/3);
h= zeros(np,1);

for j= 1:np
    uj = x(3*(j-1)+1:3*j);
    h(j) = uj'*uj-1;
end
end
