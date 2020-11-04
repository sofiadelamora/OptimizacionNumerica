% Proyecto 2 Optimización Numérica
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
figure(1) %se abre la primera figura
n=50; %número de meridiano en la esfera
np=21; %puntos


%
% Puntos aleatorios
%

figure(1)

P = randn(3,np); % matriz de 3xnp nrom. dist %P(:,j) es vector en R^3
x0 = zeros(3*np,1);  %vector que va a guardar todas las coordenadas
for j=1:np %Normalizar puntos
    v= P(:,j); % vector columna j de P (3x1)
    nv= norm(v);
    P(:,j)=v/nv; % vector normalizado
    x0(3*(j-1)+1:3*j) = P(:,j); %definimos los primeros tres elementos de x
end


% Puntos aleatorios en esfera


puntosEnEsfera(x0,"Puntos aleatorios (distribución normal)")


%
% Resultado de PCSGlobal
%

% Llamar funciones

fx = 'fesfera';     % funcion objetivo 
hx = 'hesfera';     % restricciones 


%Llamamos al método 'pcsgblobal' para encontrar el minimo


atime = cputime;
[x,lambda,k] = pcsglobal(fx,hx,x0);
timeA = cputime - atime;

figure(2)
puntosEnEsfera(x,"Mínimo obtenido utilizando la función pcsglobal.")

%
% Ahora la solución pero con Matlab
%
options.MaxFunctionEvaluations = 1.e+05;
options = optimset('Algorithm','sqp');
mtime = cputime;
[xf, fx, exitflag, output] = fmincon('fesfera',x,[],[],[],[],[],[],'hesfera_matlab');
timeM = cputime - mtime;

figure(3)
puntosEnEsfera(xf,"Mínimo obtenido en matlab")

fprintf('\n\t%s \t\t%s \t\t%s \t\t\t%s \t\t\t\t%s\n' ,'Método', 'Iteraciones','Tiempo', 'f(x)');
fprintf('\t-----------------------------------------------------------------------------------');
fprintf('\n\t%s \t\t%d \t\t\t%.4d  \t%.4d  \t\t%.4d', 'pcs_global', k, timeA, fesfera(x));
fprintf('\n\t%s \t\t%d \t\t\t\t%.4d  \t%.4d  \t\t%.4d', 'fmincon',  output.iterations, timeM,fx);
fprintf('\n');
fprintf('\t-----------------------------------------------------------------------------------');
fprintf('\n');





