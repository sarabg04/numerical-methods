function [sol, incr, incr2, iter, ACOC] = OstrowskiSist(F, dF, x0, tol, maxiter)
% OSTROWSKISIST
% Método multipaso de Ostrowski para resolver sistemas no lineales F(x)=0
% INPUT:
%   F        -> sistema no lineal F(x)=0 (función anónima)
%   dF       -> jacobiana de F(x)
%   x0       -> punto inicial
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
% Historiales de error
incr  = zeros(maxiter + 1, 1);
incr2 = zeros(maxiter + 1, 1);
incr(1)  = 1;
incr2(1) = 1;
ACOC = [];
% Forzar salida vectorial de F
F = @(x) reshape(F(x), [], 1);
% Evaluación inicial
valF  = F(x0);
valdF = dF(x0);
% Iteración principal
while incr(iter + 1) + incr2(iter + 1) > tol && iter < maxiter
    % Paso tipo Newton
    z = valdF \ valF;
    y = x0 - z;
    % Diferencia dividida
    ddF = DifDiv(F, x0, y);
    % Matriz corregida de Ostrowski
    A = 2 * ddF - valdF;
    % Evaluación en punto intermedio
    valFy = F(y);
    % Corrección final
    z = A \ valFy;
    x1 = y - z;
    iter = iter + 1;
    % Incrementos
    incr(iter + 1)  = norm(x1 - x0);
    incr2(iter + 1) = norm(F(x1));
    % Actualización
    valF  = F(x1);
    valdF = dF(x1);
    x0 = x1;
end
% Control de convergencia
if iter >= maxiter
    sol = 'No ha convergido';
else
    sol = x0;
    incr  = incr(2:iter + 1);
    incr2 = incr2(2:iter + 1);
    % Orden de convergencia
    ACOC = fACOC(incr);
end
end
