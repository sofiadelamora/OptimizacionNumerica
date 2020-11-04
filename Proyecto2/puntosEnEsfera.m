function puntosEnEsfera(cords,texto)
% cords vector con coordenadas 
% texto es una cadena de caracteres con el titulo de la gráfica

%Parametros
n = 50; %numero de meridianos en la esfera

np = floor(length(cords)/3);

[a,b,c] = sphere(n); 
w = surfl(a, b, c);       % crea la superficie con respecto a los ejes
set(w, 'FaceAlpha', 0.2);  % transparencia en la imagen
set(w, 'FaceColor', [0/255 255/255 162/255]);
axis equal   
hold on  % para añadir la gráfica

for j=1:np %Normalizar puntos
    scatter3(cords(3*j-2),cords(3*j-1),cords(3*j), 50, 'MarkerEdgeColor', [0 0 0] ,'MarkerFaceColor',[193/255 80/255 86/255]);
    hold on
end
title(texto)