function I = GaussLegendre(f,a,b,c,d,n,m)
% GAUSSLEGENDRE
% Aproximación de integrales dobles mediante cuadratura de Gauss-Legendre
% INPUT:
%   f(x,y)  -> función a integrar
%   [a,b]   -> intervalo en x
%   [c,d]   -> intervalo en y
%   n       -> número de nodos en x
%   m       -> número de nodos en y
% OUTPUT:
%   I       -> aproximación de la integral dobl
% Nodos y pesos en [-1,1] para x
[nodosx, pesosx] = NPLegendre(n);
% Cambio de intervalo [-1,1] -> [a,b]
x = (a+b)/2 + (b-a)/2 * nodosx;
% Nodos y pesos en [-1,1] para y
[nodosy, pesosy] = NPLegendre(m);
% Cambio de intervalo [-1,1] -> [c,d]
y = (c+d)/2 + (d-c)/2 * nodosy;
% Inicialización acumulador
Ix = zeros(n,1);
% Integración interna en y
for i = 1:n
    z = feval(f, x(i), y);
    Ix(i) = ((d-c)/2) * sum(pesosy .* z);
end
% Integración externa en x
I = ((b-a)/2) * sum(pesosx .* Ix);
end
