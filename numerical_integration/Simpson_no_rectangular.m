function [I]=Simpsonnorectangular(f,a,b,func,fund,n,m)
%SIMPSONNORECTANGULAR
% Aproximación de integrales dobles en dominios no rectangulares
% mediante regla compuesta de Simpson
% INPUT:
%   f     -> función f(x,y)
%   [a,b] -> intervalo en x
%   func  -> frontera inferior c(x)
%   fund  -> frontera superior d(x)
%   n     -> particiones en x
%   m     -> particiones en y
% OUTPUT:
%   I -> aproximación de la integral
% Discretización en x
h = (b - a) / (2*n);
x = a:h:b;
% Pesos de Simpson en x
pesosx = 2 * ones(1, 2*n + 1);
pesosx(1) = 1;
pesosx(2:2:end-1) = 4;
pesosx(end) = 1;
% Pesos de Simpson en y
pesosy = 2 * ones(1, 2*m + 1);
pesosy(1) = 1;
pesosy(2:2:end-1) = 4;
pesosy(end) = 1;
% Inicialización
Iy = zeros(1, 2*n + 1);
% Integración en dominio variable
for i = 1:2*n + 1
    % Fronteras dependientes de x
    c = feval(func, x(i));
    d = feval(fund, x(i));
    % Caso degenerado
    if c == d
        Iy(i) = 0;
        continue;
    end
    % Malla en y dependiente de x
    k = (d - c) / (2*m);
    y = c:k:d;
    % Evaluación de la función
    fi = feval(f, x(i), y);
    % Integración en y
    Iy(i) = k/3 * sum(pesosy .* fi);
end
% Integración final en x
I = (h/3) * sum(pesosx .* Iy);
end
