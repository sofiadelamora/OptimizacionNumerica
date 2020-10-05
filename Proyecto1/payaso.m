%Se resuelve el problema de descenso en dos pasos con diferentes valores
%para la imagen en grises del payaso
clear all;
clc;
load clown; %Carga imagen con la que se trabaja
colormap('gray');%Pone la imagen del payaso en grises
image(X); %Se abre la imagen del payaso en grises
for k = [5, 20, 30, 60, 80] %Calculamos para distintas k la factorizacion
    [W, H] = descenso2pasos(X, k);
    norm(X-W*H);
    clown=image(W*H);%Genera la imagen
    saveas(clown,sprintf('Clown_k%d.png',k));%Guarda la imagen con el nombre "Clown_k#" con la k que se utiliza
end