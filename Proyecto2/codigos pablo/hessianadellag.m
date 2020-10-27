function [ B ] = hessianadellag(fname,hname,x,y )
% Se calcula la matriz hessiana del Lagrangeano con respecto a x
%  del problema
% Min  f(x)
% s.a.  h(x) = 0
%
% donde f: R^n --> R, h: R^n --> R^m , m >0 n y f,h son de clase C^2
% 
% In
% fname .- Cadena de caracteres con la función objetivo, f(x) 
%  y la  restricción , h(x).
%          Llamado: [fx, hx] = feval(fname,x);
% x .-  Vector inicial de dimensión n.
% y.- multiplicador de Lagrange de dimensión m.
% Out
% B.- matriz nxn con la aproximación de la matriz hessiana del Lagrangeano
%     con respecto a x.
%-------------------------------------------------
%  28 febrero de 2018
% ITAM
% Optimización Numérica
%---------------------------------------------------------------

paso = 1.e-06;              % tamaño de paso para aproximar las derivadas
n = length(x);              % dimensión de la variable primal x
m = length(y);              % dimensión de las restricciones
B = zeros(n);               % matriz hessiana del lagrangiano

% [f0, h0] = feval(fname,x);  % función objetivo en el punto
f0 = feval(fname,x);
h0 = feval(hname,x);

for i = 1:n
    x1 = x;  x1(i) = x1(i) + paso;      % paso en coordenada i
%     [f1, h1] = feval(fname,x1);         % evaluación en el punto #1
    f1 = feval(fname,x1);
    h1 = feval(hname,x1);
    for j = i:n
        x2 = x; x2(j) = x2(j)+ paso;    % paso en la coordenada j
%         [f2, h2] = feval(fname,x2);     % evaluación en el punto # 2
        f2 = feval(fname,x2);
        h2 = feval(hname,x2);
        x3 = x; x3(i) = x3(i) + paso;
        x3(j) = x3(j) + paso;           % paso en coordenada i+j
%         [f3, h3] = feval(fname,x3);     % evaluación en el punto # 3
        f3 = feval(fname,x3);
        h3 = feval(hname,x3);
        B(i,j) =  (f0-f1-f2+f3);        % hessiana de f(x)
        for k =1:m                 
           B(i,j) = B(i,j) + y(k)*(h0(k)-h1(k)-h2(k)+h3(k)); 
        end
        B(i,j) = B(i,j)/(paso^2);
        if (i~=j) B(j,i) = B(i,j); end
    end
end
        

end

