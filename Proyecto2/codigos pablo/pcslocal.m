function [x,lambda,iter] = pcslocal(fname,hname,x)

% 15 de marzo de 2017

% Programacion Cuadratica Sucesiva Local
% para el problema
%    min     f(x)
%    s. a.   h(x) = 0
%
%  cuyo Lagrangeano es: L(x,lambda)= f(x) + h(x)'*lambda
%------------------------------------------------------------------
% In
% fname .- Cadena con las funciones objetivo, f(x),
%   
% hname.- Cadena con los restricciones h(x)
%          
% x .- Punto inicial
%
% 
% Out
% x .- Punto que satisface las condiciones necesarias de primer
%      orden a cierta tolerancia.
% lambda.- multiplicador de Lagrange, vector columna de dimensión m.
% iter.- número de iteraciones
%------------------------------------------------------
%             
%  * El multiplicador de Lagrange, lambda, se calcula
%    como el multiplicador de Lagrange en el subproblema cuadrático.
%  * La matriz Hessiana del Lagrangeano se calcula con el programa
%       hessianalag.m
%------------------------------------------------------
%
% Funciones en Matlab: feval.m, tril.m
% Funciones propias en Matlab:
%      gradiente.m / aproximar gradientes
%      jacobiana.m / aproxima matriz jacobiana
%      hessianalag.m / aproxima matriz hessiana del Lagrangeano
% Resultados en formato exponencial largo                     
%------------------------------------------------------
%             
% Constantes
% tol .- Tolerancia para las condiciones de primer
%        orden.
% maxiter .- Numero maximo permitido de iteraciones en PCS.
%
%------------------------------------------------------
format long e
home
tol = 1.e-06;      
maxiter = 50;

% valores iniciales
iter = 0;
hx = feval(hname,x);
df = gradiente(fname,x);
dh = jacobiana1(hname,x);
n = length(x);
m = length(hx);

% Multiplicador de Lagrange

lambda = ones(m,1);
fin = norm([[df + dh'*lambda]'  hx']);
disp(' iter     Condiciones de Primer Orden')
disp('-------------------------------------')
disp(sprintf(' %2.0f      %3.16f',iter, fin))

while ( (fin >= tol) & (iter < maxiter))
   % Hessiana con respecto a x del Lagrangeano
    B = hessianalag(fname,hname,x,lambda);
 %--------------------------------------------------
 % Solucion del Problema Cuadratico
   M = [ B  dh' ; dh zeros(m)];  %Matriz del Sistema Lineal
   ld = -[df; hx];                %Lado Derecho
   
   ps = M\ld;
   %--------------------------------------------------
   x = x + ps(1:n);      % nuevo punto
   lambda = ps(n+1:n+m); %  ¿ nuevo multiplicador de Lagrange ?
   %--------------------------------------------------------------
   hx = feval(hname,x);
   df = gradiente(fname,x);
   dh = jacobiana1(hname,x);
    fin = norm([[df + dh'*lambda]'  hx']);  
    iter = iter +1;
    disp(sprintf(' %2.0f   %3.16f',iter,fin))
end
     