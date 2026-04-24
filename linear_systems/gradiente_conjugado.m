function [sol,R,iter]=gradienteconjugado(A,b,x0,maxiter,tol)
% GRADIENTECONJUGADO
% Método del gradiente conjugado para resolver sistemas lineales Ax = b
% SOLO válido para matrices simétricas y definidas positivas
% INPUT:
%   A        -> matriz del sistema
%   b        -> vector de términos independientes
%   x0       -> aproximación inicial
%   maxiter  -> número máximo de iteraciones
%   tol      -> tolerancia de parada
% OUTPUT:
%   sol   -> solución aproximada o mensaje si no converge
%   R     -> historial del residuo ||Ax - b||
%   iter  -> número de iteraciones realizadas
iter = 0;
% Historial de residuos
R = [];
% Residuo inicial
r0 = b - A*x0;
% Dirección inicial
d0 = r0;
% Solución inicial
x1 = x0;
% Iteración principal
while iter < maxiter
    % Paso óptimo del método
    t = (r0' * r0) / (d0' * A * d0);
    % Actualización de la solución
    x1 = x1 + t * d0;
    % Nuevo residuo
    r1 = b - A * x1;
    % Norma del residuo
    residuo = norm(r1);
    % Guardar historial
    R = [R, residuo];
    % Criterio de parada
    if residuo <= tol
        sol = x1;
        return;
    end
    % Actualización del coeficiente beta
    beta = (r1' * r1) / (r0' * r0);
    % Nueva dirección de descenso
    d0 = r1 + beta * d0;
    % Actualización de variables
    r0 = r1;
    iter = iter + 1;
end
% Si no converge
if residuo > tol
    sol = 'Necesito más iteraciones :c';
else
    sol = x1;
end
end
