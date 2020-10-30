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
format long e
%Valores iniciales
n = length(x0);
m = length(hx);
k=0;
c_1=10^-02; 
Ck=1; %C_0 inicializada en 1
lambda=zeros(m,1);% Multiplicador de Lagrange
B= eye(n);%B_0 Hessiana con respecto a x del Lagrangeano
hx_eval = feval(hx,x0); %Restricciones evaluadas en la
gf = gradiente(fx,x0); %gradiente
Ak = jacobiana(hx,x0) ; %jacobiana
Cmax=10^5;
x=x0;%el vector x inicia en x_0
%----------------------------------------------------------------------
%Valores de paro
tol=10^-05;
maxk=100;
vk=[[gf + Ak'*lambda]'  hx_eval'];
%----------------------------------------------------------------------
%Metodo
while ((norm(vk) >= tol) && (k < maxk))
    
    %Resuelve problema cuadratico
    W = [ B  Ak' ; Ak zeros(m)];  %Matriz del sistema lineal
    ld = -[gf; hx_eval]; %Lado derecho
   
    sol = W\ld;
    pk = sol(1:n); %Solución en la iteración k 
    lambda = sol(n+1:n+m);%Nuevo multiplicador de Lagrange (lambda k+1)

    %------------------------------------------------------------
    %Escoger (actualizar) Ck+1
    if gf'*pk-Ck*norm(hx_eval,1)>=0
        Ck=min(Cmax, abs(gf'*pk)/norm(hx_eval,1)+1);   
    end
    
    %------------------------------------------------------------
    %Busqueda en linea
    alfak=1;
    Dk= gradiente(fx,x)' * pk - Ck* norm(feval(hx,x),1);
    max_iter=2500;
    t=0;
    while ( Fn_merito(fx,hx,x+alfak*pk,Ck) > Fn_merito(fx,hx,x,Ck)+alfak*c_1*Dk) 
       alfak=alfak/2;
       t=t+1;
    end 
        %------------------------------------------------------------
    %Guardar parámetros de la iteración k
    Ck_ant=Ck;
    x_ant=x;
    h_ant = hx_eval;
    g_ant = gf;
    A_ant = Ak;
    lambda_ant = lambda;
    %------------------------------------------------------------
    %Actualizar en k+1
    x = x + alfak*pk; %x_k+1
    s = x -x_ant; %s_k
    y = [gf + Ak'*lambda] - [gf + Ak'*lambda]; %y_k
    hx_eval = feval(hx,x); 
    %------------------------------------------------------------
    %------------------------------------------------------------
    %Actualizacion Hessiana con BFGS powell
    if s'*y <= 0.2*s'*B*s
       theta = (0.8*s'*B*s)/(s'*B*s -s'*y);
       r = y +theta*(B*s -y);
    else
       r = y;
    end
    
    B = B - (B*s*s'*B)/(s'*B*s) + (r*r')/(s'*r);
    
    if cond(B)> 10^4
        B=eye(n);
    end
   %------------------------------------------------------------
   %Nuevo multiplicador de Lagrange
    %Resolviendo min || gf_k+1 +A_k+1' * lambda||2
                %lambda en R^m
   %Solucion: lambda= -(A_k+1'*A_k+1)^-1*A_k+1*gf_k+1
   gf = gradiente(fx,x); %gradiente en k+1
   Ak = jacobiana(hx,x) ; %jacobiana en k+1
   lambda= -(Ak*Ak')\Ak*gf;
   %------------------------------------------------------------
   % Actualizar iteracion
   k=k+1;
   %------------------------------------------------------------
   % Actualizar vk
   vk=[[gf + Ak'*lambda]'  hx_eval'];
end
disp(' iteraciones     Norma de las CNPO')
disp('-------------------------------------')
disp(sprintf(' %i      %3.16f',k, norm(vk)))
end
%------------------------------------------------------------

