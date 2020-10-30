%% Proyecto 2 Optimización Numérica
%Luis Guillermo Pizana
%Sofia De la Mora 
%Optimizacion Numerica
%Proyecto 2
%05 noviembre 2020
%Se resuelve el problema 'fesfera' sujeto a las restricciones 'hesfera'
%np = 21 puntos.
%Se fija el primer punto en u1 = (1, 0, 0)'
clc;clear all;close all;

% Datos de esfera

np=21;
x = randn(3*np,1);

%Dibujos
[X, Y, Z]=sphere(50);
axis equal
hold on
%Graficar puntos
for k =1:np
    z = x(3*(k-1)+1:3*k);
    plot3(z(1), z(2), z(3),'dr','Linewidth',3)
    hold on
end
pause
close all

% Llamar funciones
fx = 'fesfera';     % funcion objetivo 
gx = 'hesfera';      % restricciones 

%Llamamos al codigo para encontrar el minimo
[x,lambda,iter] = pcsgobal(fx,x,x);
