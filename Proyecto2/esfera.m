
% Proyecto 2 Optimizaci�n Num�rica
%Luis Guillermo Pizana
%Sofia De la Mora 
%Optimizacion Numerica
%Proyecto 2
%05 noviembre 2020
%Se resuelve el problema 'fesfera' sujeto a las restricciones 'hesfera'
%np = 21 puntos.
%Se fija el primer punto en u1 = (1, 0, 0)'
clc;clear all;close all;
tic;

% Datos de esfera
figure(1) %se abre la primera figura
n=50; %n�mero de meridiano en la esfera
np=21; %puntos

%Dibujos

sphere(n);
axis equal
hold on

%Graficar puntos
P = randn(3,np); % matriz aleatoria %P(:,j) es vector en R^3
x0 = zeros(3*np,1);
for j=1:np
    v= P(:,j); 
    nv= norm(v);
    P(:,j)=v/nv;
    x(3*(j-1)+1:3*j) = P(:,j);
    %plot3(P(1,j), P(2,j), P(3,j), 'rd', 'Linewidth',3)
    hold on
end
title('Puntos en la esfera')

x0 = P(:);

% Llamar funciones
fx = 'fesfera';     % funcion objetivo 
hx = 'hesfera';     % restricciones 

%Llamamos al m�todo 'pcsgblobal' para encontrar el minimo
[x,lambda,~] = pcsglobal(fx,hx,x0);

    % Graficamos la esfera y los resultados 
figure(2)
sphere(n) %esfera unitaria 
axis equal %para que plotee efericamente 
hold on %para a�adir grafica
z = zeros(3, np);
for k =1:np
    z = x(3*(k-1)+1:3*k);
    plot3(z(1), z(2), z(3),'db','Linewidth',3)
    hold on
end
title('Minimo obtenido con pcsglobal')


% Ahora la soluci�n pero con Matlab
options.MaxFunctionEvaluations = 1.e+05;
options = optimset('Algorithm','sqp');
[xf, fx, exitflag, output] = fmincon('fesfera',x,[],[],[],[],[],[],'hesfera_matlab');

figure(3)
sphere(n) %esfera unitaria 
axis equal %para que plotee efericamente 
hold on %para a�adir grafica
y = zeros(3, np);
for k =1:np
    y = xf(3*(j-1)+1:3*j);
    plot3(y(1), y(2), y(3),'db','Linewidth',3)
    hold on
end
title('Minimo obtenido con matlab')

fx-fesfera(x)

toc;


