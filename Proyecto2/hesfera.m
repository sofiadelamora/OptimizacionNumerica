function [h] = hesfera(x)
% Funci�n de restricciones del problema de np puntos en la esfera unitaria
% de dimensi�n 3.
% Codigo clase modificado
% Luis Guillermo Piza�a
% Sof�a De la Mora
% Optimizaci�n Num�rica
% (Poyecto 2)
% Se modifico el de clase para indicar que u_1=(1,0,0).

n = length(x);
np = floor(n/3);
h = zeros(np+2,1);%Recorremos 2 para a�adir restricciones %Linea modificada despues de clase
h(1) = x(1)-1; % Agrega restriccion x(1) - 1 = 0 %Linea agregada despues de clase
h(2:3) = x(2:3);% Agrega restricciones x(2) = 0 y x(3) = 0 %Linea agregada despues de clase
for j = 2:np %Linea modificada despues de clase
    uj = x(3*(j-1)+1:3*j);
    h(j+2) = uj'*uj-1;%Recorremos 2 para a�adir restricciones %Linea modificada despues de clase
end
end
