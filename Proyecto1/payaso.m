%Se resuelve el problema de descenso en dos pasos con diferentes valores
%para la imagen en grises del payaso
clear all;
clc;
load clown; %Carga imagen con la que se trabaja
colormap('gray');%Pone la imagen del payaso en grises
image(X);
%Para descenso2pasos. Descomentar lineas 37 y 48 de descenso2pasos y
%comentar las lineas 38 y 49
%Generamos las imagenes y las guardamos:
for k = 80 %[5, 20, 30, 60, 80] %Calculamos para distintas k la factorizacion
    tic;
    [W, H] = descenso2pasos(X, k);
    norm(X-W*H);
    clown=image(W*H);%Genera la imagen
    saveas(clown,sprintf('Clown2p_k%d.png',k));%Guarda la imagen con el nombre "Clown2p_k#" con la k que se utiliza
    toc;
end
%Para Quadprog. Descomentar lineas 38 y 49 de descenso2pasos y
%comentar las lineas 37 y 48
%Generamos las imagenes y las guardamos:
% for k = [5, 20, 30, 60, 80] %Calculamos para distintas k la factorizacion
%     tic;
%     [W, H] = descenso2pasos(X, k);
%     norm(X-W*H);
%     clown=image(W*H);%Genera la imagen
%     saveas(clown,sprintf('ClownQuad_k%d.png',k));%Guarda la imagen con el nombre "ClownQuad_k#" con la k que se utiliza
%     toc;
%end