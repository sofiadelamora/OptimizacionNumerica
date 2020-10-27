function [x, lambda, k] = pcsglobal(fx, hx, x0)
% Metodo de programacion Cuadratica Sucesiva con busqueda de lineal,
% usando la funcion de merito L-1 y actualizacion de la hessiana
% con la formula BF GS para el problema
% Min fx
% Sujeto a hx = 0
%
% fx y hx son cadenas de caracteres con las funciones en Matlab
% de la funcion objetivo y las restricciones del problema
% El vector x0 es el valor inicial
% Salida
% x.- aproximacion al minimo local
% lambda.- multiplicador de Lagrange asociado a x.
% k.- numero de iteraciones realizadas.
%
% Debe usar las funciones: gradiente.m y jacobiana.m para calcular
% las primeras derivadas.
%-----------------------------------------------------
%Luis Guillermo Pizana
%Sofia De la Mora 
%Optimizacion Numerica
%Proyecto 2
%05 noviembre 2020
%-----------------------------------------------------------------------
%Valores iniciales
n = length(x0);
m = length(hx);
k=0;
c_1=10e^-02;
C0=1;
lambda=zeros(m,1);% Multiplicador de Lagrange
B_0= eye(n);% Hessiana con respecto a x del Lagrangeano
x = zeros(n,1);
hx = feval(hx,x0);
gf = gradiente(fx,x0); %gradiente
Ak = jacobiana(hx,x0) ; %jacobiana
alfak=1;
%----------------------------------------------------------------------
%Valores de paro
tol=10e^-05;
maxk=100;
r=norm([[gf + Ak'*lambda]'  hx']);
%----------------------------------------------------------------------
%Metodo
while ((fin >= tol) & (iter < maxiter) & (deltak >= tol)
    %Resuelve problema cuadratico
    W = [ B  Ak' ; Ak zeros(m)];  %Matriz del Sistema Lineal
    ld = -[gf; hx]; %Lado Derecho
   
    sol = W\ld;
    pk = sol(1:n);
    
    %Escoger Ck+1
    if gf'*pk/norm(hx,1)>0
        Ck= gf'*pk/norm(hx,1);
    elseif gf'*pk/norm(hx,1)== 0
        Ck=1;
    else 
        Ck=-gf'*pk/norm(hx,1);
    end
        
    %Iteracion anteior
    x_ant = x;
    h_ant = hx;
    g_ant = gf;
    A_ant = dh;
    lambda_ant = lambda;
    
    %Actualizar x
    x = x+alfak*pk;
    lambda = ps(n+1:n+m);  %Nuevo multiplicador de Lagrange (lambda k+1)
    
    hx = feval(hx,x);
    gf = gradiente(fx,x);
    Ak = jacobiana1(hx,x);
    
    %Recorte de paso
    phik = feval(fx,x_ant)+feval(hx,x_ant)'*lambda+Ck/2*norm(feval(hx,x_ant))^2;
    phik_uno = feval(fx,x)+feval(hx,x)'*lambda+Ck/2*norm(feval(hx,x))^2;
    Dk= feval(gf, x_ant)' * pk - Ck* norm(feval(hx,x),1);
    while (phik_uno> phik + alfak * Dk)
       alfak=alfak/2; 
    end
    
    s = x -x_ant;
    y = [g_ant + A_ant'*lambda] - [gf + Ak'*lambda];
    %------------------------------------------------------------
    %Actualizacion Hessiana con BFGS powell
    if s'*y <= 0.2*s'*B*s
       theta = (0.8*s'*B*s)/(s'*B*s -s'*y);
       r = y +theta*(B*s -y);
    else
       r = y;
    end
    
    if cond(B)> 10^4
        B=eye(n);
    else
        B = B - (B*s*s'*B)/(s'*B*s) + (r*r')/(s'*r);
    end
    
end
end