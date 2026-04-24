function [sol,iter,incr,incr2,ACOC]=newton(F,DF,x0,tol,maxiter)
% NEWTON
% Método de Newton para sistemas no lineales F(x)=0
% INPUT:
%   F        -> función vectorial F(x)
%   DF       -> jacobiana de F(x)
%   x0       -> punto inicial
%   tol      -> tolerancia
%   maxiter  -> número máximo de iteraciones
% OUTPUT:
%   sol    -> solución aproximada o mensaje de no convergencia
%   iter   -> número de iteraciones
%   incr   -> historial de incrementos ||x_k - x_{k-1}||
%   incr2  -> historial de residuos ||F(x_k)||
%   ACOC   -> orden de convergencia estimado
% Asegurar vector columna
x0 = x0(:);
incr  = tol + 1;
incr2 = tol + 1;
iter = 0;
% Historial de incrementos
I = [];
% Iteración de Newton
while incr + incr2 > tol && iter < maxiter
    % Evaluación de función y jacobiana
    Fx = F(x0);
    Jx = DF(x0);
    % Paso de Newton
    z = Jx \ (-Fx);
    % Actualización
    x1 = x0 + z;
    % Errores
    incr  = norm(x1 - x0);
    incr2 = norm(F(x1));
    % Guardar historial
    I = [I incr];
    % Preparar siguiente iteración
    x0 = x1;
    iter = iter + 1;
end
% Control de convergencia
if incr > tol
    sol = 'se necesitan más iteraciones :c';
else
    sol = x0;
end
% Orden de convergencia
ACOC = log(I(3:end)./I(2:end-1)) ./ log(I(2:end-1)./I(1:end-2));
end
