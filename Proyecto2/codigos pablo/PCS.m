function [x, iter] = PCS(fname, hname, x)
% Programacion Cuadratica Sucesiva usando el
% para el problema
% min f(x)
% s. a. h(x) = 0
% In
%   fname.- Cadena de caracteres con la funciónn objetivo, f(x).
%   hname.- Cadena de caracteres con la función de restricciones h(x) = 0
%   h : Rn --> Rm.
%   x.- Punto inicial
%
% Out
%   x.- Punto que satisface las condiciones necesarias de primer
%   orden a cierta tolerancia.
%   iter.- número de iteraciones usadas en el método
%---------------------------------------------------------------------
%% Declaraciones iniciales
tol = 1.e-05;
maxiter = 100;
iter = 0;

%% Paso 1 Escoger (x0,lambda0) usando Newton
gk = gradiente(fname,x); %nos devuelve una función
hk = feval(hname, x);
Ak = jacobiana1(hname,x); %nos devuelve una función
n = length(hk);
lambda = ones(n,1);
v1 = gk + Ak'*lambda;
v2 = hk ;
F = [v1; v2];

%% Paso 2 Resolver Problema y actualizar
    while ( tol < norm(F,2) && iter < maxiter)
        Bk  = hessianadellag(fname,hname,x,lambda);
        % Resolvemos el problema de minimizar
        [pk, lambdak] = pcdirecto(Bk,Ak,gk,-hk);
        x = x + pk;
        lambda = lambdak;
        iter = iter + 1;
        
        %Actulizamos
        gk = gradiente(fname,x); %nos devuelve una función
        hk = feval(hname, x);
        Ak = jacobiana1(hname,x); %nos devuelve una función
        v1 = gk + Ak'*lambda;
        v2 = hk ;
        F = [v1; v2];
    end
end