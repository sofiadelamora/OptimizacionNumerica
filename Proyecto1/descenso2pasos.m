function [W, H] = descenso2pasos(X, k)
%Metodo de descenso en dos pasos o descenso por coordenadas
% Min   ||X - W*H ||_F^2
% s.a.   W >= 0
%        H >= 0
%
% Input
% X matriz real de rxp con entradas positivas
% k iteraciones del metodo
%
%Output
% W matriz de rxk con entradas positivas
% H matriz de kxp con entradas positivas
% W*H aproxima la matriz X
%
%--------------------------------------------------------------------------
% Optimización Numérica
% Proyecto 1
% 08 de octubre del 2020
% Luis Guillermo Pizana
% Sofía De la Mora
%--------------------------------------------------------------------------
% Parametros iniciales
[r,p] = size(X);    
W = ones(r,k);      
H = one(k,p);      




