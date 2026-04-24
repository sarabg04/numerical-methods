function [sol,iter,Acoc,GS]=GaussSeidel(A,b,x0,tol,maxiter)
% GAUSSSEIDEL
% Método iterativo de Gauss-Seidel para resolver sistemas lineales Ax=b
% INPUT:
%   A        -> matriz del sistema
%   b        -> vector de términos independientes
%   x0       -> aproximación inicial
%   tol      -> tolerancia
%   maxiter  -> número máximo de iteraciones
% OUTPUT:
%   sol   -> solución aproximada o mensaje si no converge
%   iter  -> número de iteraciones realizadas
%   Acoc  -> estimación de orden de convergencia
%   GS    -> matriz de iteración del métod
iter = 0;
incr = tol + 1;
% Descomposición de A
L = tril(A,-1);   % parte inferior estricta
U = triu(A,1);    % parte superior estricta
% Diagonal
invD = diag(1./diag(A));
D = inv(invD);
% Matriz auxiliar
matsum = D + L;

% Matriz de iteración Gauss-Seidel
GS = -inv(D + L) * U;
% Historial de errores
Ac = [];
while iter <= maxiter && incr > tol
    x = matsum \ (b - U*x0);
    Ac = [Ac, norm(x - x0)];
    incr = norm(x - x0);
    iter = iter + 1;
    x0 = x;
end
% Control de convergencia
if incr > tol
    sol = 'necesito más iteraciones :c';
else
    sol = x;
end
% Estimación del orden de convergencia
Acoc = (log(Ac(3:end)./Ac(2:end-1))) ./ (log(Ac(2:end-1)./Ac(1:end-2)));
% Interpretación:
% p ≈ 1 → convergencia lineal
% p ≈ 2 → convergencia cuadrática, etc
end
