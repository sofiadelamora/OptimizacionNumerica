function [x, lambda, k] = pcsglobal(fx, hx, x0)
% M´etodo de programaci´on Cuadr´atica Sucesiva con b´usqueda de lineal,
% usando la funci´on de m´erito L-1 y actualizaci´on de la hessiana
% con la f´ormula BF GS para el problema
% Min fx
% Sujeto a hx = 0
%
% fx y hx son cadenas de caracteres con las funciones en Matlab
% de la funci´on objetivo y las restricciones del problema
% El vector x0 es el valor inicial
% Salida
% x.- aproximaci´on al m´?nimo local
% ?.- multiplicador de Lagrange asociado a x.
% k.- n´umero de iteraciones realizadas.
%
% Debe usar las funciones: gradiente.m y jacobiana.m para calcular
% las primeras derivadas.
%
% Nombres completos de los integrantes del equipo
% Instrucciones del programa documentados.
end