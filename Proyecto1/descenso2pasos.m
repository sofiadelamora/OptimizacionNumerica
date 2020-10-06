function [W, H] = descenso2pasos(X, k)
%Metodo de descenso en dos pasos o descenso por coordenadas para
%factorización no negativa de una matriz
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
% Optimización Numérica
% Proyecto 1
% 08 de octubre del 2020
% Luis Guillermo Pizana
% Sofía De la Mora
%--------------------------------------------------------------------------
% Parametros iniciales
[r,p] = size(X);    
W = ones(r,k);      
H = ones(k,p);
%--------------------------------------------------------------------------
%
for i = 1:5%Comienzan las iteraciones. Se fijaron el máximo de iteraciones 
    %en 5 porque a partir de la sexta los cambios son mínimos.
    %Calcula H fijando W
    for j = 1:p
        Qh= W'*W;
        ch = (-X(:, j)'*W)';
        bh = zeros(k,1);
        Ah = eye(k);
        [xh,~,~] = punintpc(Qh, Ah, ch, bh);
        %xh=quadprog(Qh,ch,-Ah,-bh);
        H(:,j) = xh; 

    end
    % Calcula W fijando H
    for l = 1:r
        Qw= H*H';
        cw = (-X(l, :)*H')'; 
        bw = zeros(k,1);
        Aw = eye(k);
        [xw,~,~] = punintpc(Qw, Aw, cw, bw);
        %xw=quadprog(Qw,cw,-Aw,-bw);
        W(l,:)= xw';

    end
end 

