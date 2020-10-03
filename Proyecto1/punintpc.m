function [x, y, mu] = punintpc(Q, A, c, b)
% Resuelve el problema cuadratico por el metodo de punto interior
% Min   (1/2)* x' * Q * x + c'* x
% s.a.   A * x >= b
%
% Input:
%
% Q matriz de nxn simétrica y definida positiva
% A matriz de mxn de rango completo con n>=m
% c es vector columna real de tamaño n
% b es vector columna real de tamaño m
%
%Output:
%
% x vector real de tamaño n que aproxima la solucion del minimo local del
%problema
% y vector de tamaño m que contiene las variables de holgura del problema
% mu vector de tamaño m que contiene el multiplicador de lagrange de la
%restricción desigual
%
%Supuesto:
%
% Suponemos que el conjunto de vectores x en R^n con Ax>b es distinto al
% vacío
%
%--------------------------------------------------------------------------
% Optimización Numérica
% Proyecto 1
% 08 de octubre del 2020
% Luis Guillermo Pizana
% Sofía De la Mora
%--------------------------------------------------------------------------
% Parametros iniciales
n = length(c);     % numero de variables a minimizar
m = length(b);     % numero de restricciones
maxiter = 250;     % iteraciones maximas
tol = 1e-08;      % tolerancia CNPO 
citer = 0;         % contador iteraciones
%--------------------------------------------------------------------------
%Variables iniciales
x = ones(n,1); 
y = ones(m,1);
mu = ones(m,1);
e = ones(m,1);
eta = (0.5)*(y'*mu)/m; %sigma*y'mu/m, con sigma en [0,1] tomamos 1/2.
%--------------------------------------------------------------------------
H =[Q*x - A'*mu+c; A*x- y - b; mu.*y];
norma = norm(H); %Norma de CNPO
while(norma > tol && citer < maxiter)
    % Resuelve el sistema lineal de Newton para la trayectoria central
    Y=diag(y);
    U=diag(mu);
 
    r_x=Q*x-A'*mu+c; 
    r_y=A*x-y-b;     
    r_mu=Y*U*e-eta*e; 

    %Sistema 
    D= Q+A'*diag(mu./y)*A;
    z=-(r_x+A'*diag(mu./y)*r_y+A'*(r_mu./y));
    
    %Pasos
    delta_x= z\D;
    delta_y=A*delta_x'+r_y;
    delta_mu=-(diag(mu./y)*delta_y-r_mu./y);
    
    %Recorte de paso
    bt = []; gm = [];
    for k =1:m
        if(delta_mu(k) < 0)
            bt = [bt; -(mu(k)/delta_mu(k))];
        else
            bt= [bt; 1];
        end
        if (delta_y(k) < 0)
            gm = [gm; -(y(k)/delta_y(k))];
        else
            gm = [gm; 1];
        end
    end
    
    alfa = min([bt ; gm]);
    alfa =(0.9995)*min([1 alfa]);  
    
    % Nuevo punto
      
    x = x + alfa*delta_x';
    mu = mu + alfa*delta_mu;
    y = y + alfa*delta_y;
    %Nueva eta
    eta = (0.5)*(y'*mu)/m;
%     if(mod(citer,2)==0)
%        eta = (0.5)*(y'*mu)/m;
%      else
%        eta = 0;
%      end
    %CNPO
    H =[Q*x - A'*mu+c; A*x - y - b; mu.*y];    
    norma = norm(H);
    citer = citer + 1;
end
end

