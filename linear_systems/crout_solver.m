function sol = crout_solver(a, b, c, d)
% CROUT_SOLVER
% Resuelve un sistema tridiagonal Ax = d usando factorización LU (Crout)
% INPUT:
%   a - diagonal principal
%   b - diagonal superior
%   c - diagonal inferior
%   d - vector independiente
% OUTPUT:
%   sol - solución del sistema Ax = d
n = length(a);
% Inicialización
    l = zeros(n,1);
    u = zeros(n-1,1);
    z = zeros(n,1);
    x = zeros(n,1);
% Factorización LU (Crout)
    l(1) = a(1);
    u(1) = b(1) / l(1);
    for i = 2:n-1
        l(i) = a(i) - c(i-1) * u(i-1);
        u(i) = b(i) / l(i);
    end
    l(n) = a(n) - c(n-1) * u(n-1);
% Resolución Lz = d (sustitución hacia adelante)
    z(1) = d(1) / l(1);
    for i = 2:n
        z(i) = (d(i) - c(i-1) * z(i-1)) / l(i);
    end
% Resolución Ux = z (sustitución hacia atrás)
    x(n) = z(n);
    for i = n-1:-1:1
        x(i) = z(i) - u(i) * x(i+1);
    end
    sol = x;
end
