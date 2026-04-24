function [sol,iter,R,Acoc,J] = Jacobi(A,b,x0,maxiter,tol)
% JACOBI
% Método iterativo de Jacobi para resolver sistemas lineales Ax = b
% INPUT:
%   A        -> matriz del sistema
%   b        -> vector de términos independientes
%   x0       -> aproximación inicial (vector)
%   maxiter  -> número máximo de iteraciones
%   tol      -> tolerancia
% OUTPUT:
%   sol   -> solución aproximada o mensaje de no convergencia
%   iter  -> número de iteraciones
%   R     -> historial del residuo ||Ax - b||
%   Acoc  -> estimación del orden de convergencia
%   J     -> matriz de iteración de Jacobi
iter = 0;
residuo = tol + 1;
% Historiales
R  = [];  % residuos
Ac = [];  % incrementos
% Descomposición de A
L = tril(A,-1);   % parte inferior estricta
U = triu(A,1);    % parte superior estricta
% Inversa de la diagonal (más eficiente que inv)
invD = diag(1./diag(A));
% Matriz de iteración Jacobi
J = -invD * (L + U);
% Vector constante
C = invD * b;
% Iteración principal
while iter < maxiter && residuo > tol
    % Paso de Jacobi
    x = J * x0 + C;
    % Incremento entre iteraciones
    Ac = [Ac, norm(x - x0)];
    iter = iter + 1;
    % Residuo del sistema
    residuo = norm(A*x - b);
    % Guardar historial
    R = [R, residuo];
    % Actualizar solución
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
