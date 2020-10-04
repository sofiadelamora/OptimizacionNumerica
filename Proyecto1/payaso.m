clear all;
clc;
load clown;
colormap('gray');
image(X);
for k = [5, 20, 30, 60, 80]
    [W, H] = descenso2pasos(X, k);
    norm(X-W*H)
end