function [sol, incr, incr2, iter, ACOC] = jarratt(F, dF, x0, tol, maxiter)
% JARRATT
% Método multipaso de Jarratt para resolver sistemas no lineales F(x)=0
% INPUT:
%   F        -> sistema no lineal (F(x)=0), función anónima
%   dF       -> jacobiana de F(x)
%   x0       -> aproximación inicial
%   tol      -> tolerancia
%   maxiter  -> número máximo de iteraciones
% OUTPUT:
%   sol    -> solución aproximada o mensaje de no convergencia
%   incr   -> historial de incrementos ||x_k - x_{k-1}||
%   incr2  -> historial de residuos ||F(x_k)||
%   iter   -> número de iteraciones
%   ACOC   -> orden de convergencia estimado
% Asegurar vector columna
x0 = x0(:);
iter = 0;
% Historiales de convergencia
incr  = zeros(maxiter + 1, 1);
incr2 = zeros(maxiter + 1, 1);
incr(1) = 1;
incr2(1) = 1;
ACOC = [];
% Forzar salida vectorial
F = @(x) reshape(F(x), [], 1);
% Evaluación inicial
valF  = F(x0);
valdF = dF(x0);
% Iteración principal
while incr(iter + 1) + incr2(iter + 1) > tol && iter < maxiter
    % Paso 1 del método
    z = (2/3) * (valdF \ valF);
    y = x0 - z;
    % Paso 2
    z = valdF \ valF;
    valdFy = 3 * dF(y);
    c = (valdFy + valdF) * z;
    % Corrección Jarratt
    z = (1/2) * ((valdFy - valdF) \ c);
    % Nueva aproximación
    x1 = x0 - z;
    iter = iter + 1;
    % Incremento entre iteraciones
    incr(iter + 1) = norm(x1 - x0);
    % Actualización de funciones
    valF  = F(x1);
    valdF = dF(x1);
    % Residuo
    incr2(iter + 1) = norm(valF);
    x0 = x1;
end
% Control de convergencia
if iter >= maxiter
    sol = 'No ha convergido';
else
    sol = x0;
    incr  = incr(2:iter + 1);
    incr2 = incr2(2:iter + 1);
    % Orden de convergencia estimado
    ACOC = fACOC(incr);
end
end
