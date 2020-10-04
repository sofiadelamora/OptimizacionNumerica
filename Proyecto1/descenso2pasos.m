function [W, H] = descenso2pasos(X, k)
%Metodo de descenso en dos pasos o descenso por coordenadas para
%factorizaci�n no negativa de una matriz
% Min   ||X - W*H ||_F^2
% s.a.   W >= 0
%        H >= 0
%
% Input
% X matriz real de rxp con entradas positivas
% k iteraciones del metodo, es mucho menor que r y que p
%
%Output
% W matriz de rxk con entradas positivas
% H matriz de kxp con entradas positivas
% W*H aproxima la matriz X
%
%--------------------------------------------------------------------------
% Optimizaci�n Num�rica
% Proyecto 1
% 08 de octubre del 2020
% Luis Guillermo Pizana
% Sof�a De la Mora
%--------------------------------------------------------------------------
% Parametros iniciales
[r,p] = size(X);    
W = ones(r,k);      
H = ones(k,p);
%--------------------------------------------------------------------------
%Calcula H y W
for i = 1:k
    %Calcula H fijando W
    for j = 1:p
        Qh= W'*W;
        ch = -X(:, j)'*W;
        bh = zeros(k,1);
        Ah = eye(k);
        [xh,~,~] = punintpc(Qh, Ah, ch', bh);
        H(:,j) = xh; %Creo que aqu� esta el error 1

    end
    % Calcula W fijando H
    for l = 1:r
        Qw= H'*H;
        cw = -X(l, :)*H'; 
        bw = zeros(k,1);
        Aw = eye(k);
        [xw,~,~] = punintpc(Qw, Aw, cw', bw);
        W(l,:)= xw';%Creo que aqu� esta el error 2

    end
end 

