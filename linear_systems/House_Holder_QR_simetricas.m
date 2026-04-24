function [Q,R] = House_Holder_QR_simetricas(A)
% HOUSE_HOLDER_QR_SIMETRICAS
% Factorización QR mediante reflexiones de Householder
% optimizada para matrices simétricas
% INPUT:
%   A -> matriz a factorizar
% OUTPUT:
%   Q -> matriz ortogonal
%   R -> matriz triangular superior
[m,n] = size(A);
% Matrices auxiliares
vp = zeros(m,n);
v  = zeros(m,n);
% Inicialización de Q
Q = eye(m,n);
% Iteración de Householder
for k = 1:n-2
    % Construcción del vector de reflexión
    vp(1,k) = A(k+1,k) + norm(A(k+1:n,k),2);
    vp(2:n-k,k) = A(k+2:n,k);
    v(k+1:n,k) = vp(1:n-k,k);
    % Matriz de Householder
    H = eye(n) - 2*(v(:,k)*v(:,k)')/(norm(v(:,k),2)^2);
    % Caso degenerado (columna ya anulada)
    if sum(A(k:n,k)) == 0
        H = eye(n);
    end
    % Transformación simétrica
    A = H * A * H;
    % Acumulación de transformaciones ortogonales
    Q = Q * H;
end
% Matriz triangular superior resultante
R = A;
end
