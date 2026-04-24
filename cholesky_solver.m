function [L, x] = cholesky_solver(A, b)
% CHOLESKY_SOLVER Resuelve Ax = b usando factorización de Cholesky
% Requiere A simétrica y definida positiva

    % Validaciones básicas
    if ~isequal(A, A')
        error('La matriz A no es simétrica');
    end

    n = length(b);
    L = zeros(n, n);

    % Factorización de Cholesky
    for i = 1:n
        for j = 1:i
            if i == j
                L(i,j) = sqrt(A(i,i) - sum(L(i,1:j-1).^2));
            else
                L(i,j) = (A(i,j) - sum(L(i,1:j-1).*L(j,1:j-1))) / L(j,j);
            end
        end
    end

    % Sustitución hacia adelante (Ly = b)
    y = zeros(n,1);
    for i = 1:n
        y(i) = (b(i) - sum(L(i,1:i-1).*y(1:i-1))) / L(i,i);
    end

    % Sustitución hacia atrás (L'x = y)
    x = zeros(n,1);
    for i = n:-1:1
        x(i) = (y(i) - sum(L(i+1:n,i).*x(i+1:n))) / L(i,i);
    end
end