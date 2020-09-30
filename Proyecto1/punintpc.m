function [x, y, mu] = punintpc(Q, A, c, b)
% Resuelve el problema cuadratico por el metodo de punto interior
% Min   (1/2)* x' * Q * x + c'* x
% s.a.   A * x >= b
% Q matriz de nxn sim�trica y definida positiva
% A matriz de mxn 
% c es un vector real de tama�o n
% Suponemos que el conjunto de vectores x en R^n con Ax>b es distinto al
% vac�o
%


