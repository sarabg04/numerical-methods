function y = diferencias_finitas(p, q, r, a, b, alfa, beta, N)
% DIFERENCIAS_FINITAS
% Resuelve el problema de contorno:
%   y'' = p(x)y' + q(x)y + r(x),  x ∈ [a, b]
%   y(a) = alfa, y(b) = beta
% usando diferencias finitas.
% INPUT:
%   p, q, r - funciones
%   a, b    - intervalo
%   alfa    - condición en x = a
%   beta    - condición en x = b
%   N       - número de puntos interiores
% OUTPUT:
%   y - solución aproximada en la malla (incluye extremos)
    h = (b - a) / (N + 1);
% Malla (incluyendo extremos)
    x = linspace(a, b, N + 2);
    x = x(2:end-1); % puntos interiores
% Inicialización
    dp = zeros(N, 1);     % diagonal principal
    ds = zeros(N-1, 1);    % diagonales secundaria
    d  = zeros(N, 1);      % término independiente
    for i = 1:N
        xi = a + i*h;
        dp(i) = 2 + h^2 * q(xi);
        if i > 1
            ds(i-1) = -1 + (h/2) * p(xi);
        end
        if i == 1
            d(i) = -h^2 * r(xi) + (1 + (h/2) * p(xi)) * alfa;
        elseif i == N
            d(i) = -h^2 * r(xi) + (1 - (h/2) * p(xi)) * beta;
        else
            d(i) = -h^2 * r(xi);
        end
    end
% Sistema tridiagonal
    A = spdiags([ds dp ds], -1:1, N, N);
    y_interior = A \ d;
% Añadir condiciones de frontera
    y = [alfa; y_interior; beta];
end
