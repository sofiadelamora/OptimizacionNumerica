function [x, lambda] = pcdirecto(Q,A,c,b)
% Solucion del problema
% Min(1/2) *x'*Q*x + c'x
% s.a. A*x = b
% Q es simetrica postiva definida

% A es mxm con m<n y de rango m
%-----------------------------------------

m =  length(b);
n = length(c);
K = [Q A'; A zeros(m)];

z =  K \ [-c;b];

x = z(1:n);
lambda = z(n+1:n+m);
end