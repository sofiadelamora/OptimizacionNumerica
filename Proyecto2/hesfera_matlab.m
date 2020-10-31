function [g,h] = hesfera(x)
% Función de restricciones del problema de np puntos en la esfera unitaria
% de dimensión tres.
% Codigo de clase
% Luis Guillermo Pizaña
% Sofía De la Mora
% Optimización Numérica
% ITAM
% 20 de octubre de 2020

g = []; % (uso de fmincon.m en Matlab)

n = length(x);
np = floor(n/3);
h = zeros(np+2,1);%Recorremos 2 para añadir restricciones
h(1) = x(1)-1; % Agrega restriccion x(1) - 1 = 0
h(2:3) = x(2:3);% Agrega restricciones x(2) = 0 y x(3) = 0
for j = 2:np
    uj = x(3*(j-1)+1:3*j);
    h(j+2) = uj'*uj-1;%Recorremos 2 para añadir restricciones
end
end
