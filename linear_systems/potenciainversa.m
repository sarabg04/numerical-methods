function [linv,z,iter,incr]=potenciainversa(A,z,tol,maxiter)
% POTENCIAINVERSA
% Método de la potencia inversa para aproximar autovalores y autovectores
% INPUT:
%   A        -> matriz cuadrada
%   z        -> vector inicial
%   tol      -> tolerancia
%   maxiter  -> número máximo de iteraciones
% OUTPUT:
%   linv  -> aproximación del autovalor dominante de A^{-1}
%   z     -> autovector asociado
%   iter  -> número de iteraciones
%   incr  -> último increment
incr = tol + 1;
iter = 0;
% Asegurar vector columna
z = z(:);
% Normalización inicial
z = z / norm(z);
% Resolución inicial del sistema
p = A \ z;
% Iteración principal
while incr > tol && iter < maxiter
    % Aproximación del autovalor de A^{-1}
    linv = z' * p;
    % Normalización del vector
    z = p / norm(p);
    % Resolver sistema lineal
    p = A \ z;
    % Criterio de parada
    incr = norm(z - p / norm(p));
    iter = iter + 1;
end
% Inverso para obtener autovalor de A
linv = linv^(-1);
% Mensaje de convergencia
if incr > tol
    disp("se necesitan más iteraciones");
end
end
