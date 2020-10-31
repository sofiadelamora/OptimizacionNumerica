function [g,h] = hesfera(x)
% Funci�n de restricciones del problema de np puntos en la esfera unitaria
% de dimensi�n tres.
% Codigo de clase
% Luis Guillermo Piza�a
% Sof�a De la Mora
% Optimizaci�n Num�rica
% ITAM
% 20 de octubre de 2020

g = []; % (uso de fmincon.m en Matlab)

n = length(x);
np = floor(n/3);
h = zeros(np+2,1);%Recorremos 2 para a�adir restricciones
h(1) = x(1)-1; % Agrega restriccion x(1) - 1 = 0
h(2:3) = x(2:3);% Agrega restricciones x(2) = 0 y x(3) = 0
for j = 2:np
    uj = x(3*(j-1)+1:3*j);
    h(j+2) = uj'*uj-1;%Recorremos 2 para a�adir restricciones
end
end
