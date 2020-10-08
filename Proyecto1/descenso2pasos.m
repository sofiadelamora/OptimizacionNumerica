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
W = zeros(r,k)+1/2;      
H = ones(k,p);
maxiter=k;
citer=0;
W_ant = zeros(r,k);%matrix auxiliar para guardar la iteracion anterior     
H_ant = zeros(k,p);%matrix auxiliar para guardar la iteracion anterior
tol=1e-08; %tolerancia para las condiciones de paro
norma = norm(W - W_ant,'fro') + norm(H - H_ant,'fro');%para condici�n de paro
%--------------------------------------------------------------------------
while(norma > tol && citer < maxiter)
    %Calcula H fijando W
    H_ant = H;
    W_ant = W;
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
    norma = norm(W - W_ant,'fro') + norm(H - H_ant,'fro');
    citer=citer+1;
    
end 
fprintf("\nIteraciones: %i\n", citer);
