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
%------------------------------------------------------
%Parametros iniciales
%n
%m
tol=10e^-05;
maxk=100;
c_1=10e^-02;
c_0=1;
lambda=zeros(m,1);
B_0= eye(n);
k=0;
end