function [phi] = merit_fn(fx, hx, x, Ck)
%funcion de merito L-1 
phi = feval(fx, x) + Ck*norm(feval(hx,x),1);
end