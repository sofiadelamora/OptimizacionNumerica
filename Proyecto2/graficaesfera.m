% graficaesfera.m
% 20 de octubre de 2020
close all;

sphere(50) %n es el número de meridiano en la esfera
axis equal %para que la esfera no esté achatada
hold on

np=30;

P=randn(3,np); % una matriz aleatoria
               % P(:,j) es un vector en R^3
x=zeros(3*np,1);
for j=1:np
    v= P(:,j); nv= norm(v);
    P(:,j)=v/nv;
    x(3*(j-1)+1:3*j) = P(:,j);
    plot3(P(1,j), P(2,j), P(3,j), 'rd', 'Linewidth',3)
    hold on
end

% tic; g=gradiente('fesfera',x); toc
% h= feval('hesfera',x);
% pause
% close all
% stem(h)

% Solución por Matlab

[xf, fx, exitflag, output] = fmincon('fesfera',x,[],[],[],[],[],[],'hesfera');
