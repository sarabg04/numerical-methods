function [sol,R,Acoc,iter,Ac]=gradiente(A,b,x0,maxiter,tol)
% GRADIENTE
% Método del gradiente para resolver sistemas lineales Ax = b
% SOLO válido para matrices simétricas y definidas positivas
% INPUT:
%   A        -> matriz del sistema
%   b        -> vector de términos independientes
%   x0       -> aproximación inicial
%   maxiter  -> número máximo de iteraciones
%   tol      -> tolerancia
% OUTPUT:
%   sol   -> solución aproximada o mensaje de no convergencia
%   R     -> historial del residuo ||Ax - b||
%   Acoc  -> estimación del orden de convergencia
%   iter  -> número de iteraciones
%   Ac    -> historial de incrementos ||x_k - x_{k-1}||
iter = 0;
residuo = tol + 1;
% Historiales
R = [];   % residuos
Ac = [];  % incrementos
% Iteración principal
while iter < maxiter && residuo > tol
    % Residuo actual
    r = b - A*x0;
    % Paso óptimo del gradiente
    t = (r' * r) / (r' * A * r);
    % Actualización
    x = x0 + t * r;
    % Incremento entre iteraciones
    Ac = [Ac, norm(x - x0)];
    iter = iter + 1;
    % Residuo euclídeo
    residuo = norm(A*x - b);
    % Guardar historial del residuo
    R = [R, residuo];
    % Actualizar punto
    x0 = x;
end
% Control de convergencia
if residuo > tol
    sol = 'necesito más iteraciones :c';
else
    sol = x;
end
% Orden de convergencia estimado
Acoc = (log(Ac(3:end)./Ac(2:end-1))) ./ (log(Ac(2:end-1)./Ac(1:end-2)));
% Interpretación:
% ≈1 → convergencia lineal
% ≈2 → cuadrática
end
