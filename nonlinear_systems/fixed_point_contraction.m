function results = fixed_point_contraction(fx, fy, n)
% FIXED_POINT_CONTRACTION
% Realiza iteraciones de punto fijo para un sistema en R^2
% y analiza si parece ser una contracción.
% INPUT:
%   fx - función que define x_{k+1} = fx(y_k)
%   fy - función que define y_{k+1} = fy(x_k)
%   n  - número de iteraciones
% OUTPUT:
%   results - matriz (n x 2) con las iteraciones [x_k, y_k]
% Valores iniciales
    x = 1;
    y = 0;
% Almacenar resultados
    results = zeros(n, 2);
% Iteraciones de punto fijo
    for i = 1:n
        x_new = fx(y);
        y_new = fy(x);
        results(i, :) = [x_new, y_new];
        x = x_new;
        y = y_new;
    end
% Mostrar resultados
    disp('Iteraciones:');
    disp(results);
% Calcular distancias entre iteraciones consecutivas
    distances = vecnorm(diff(results), 2, 2);
    disp('Distancias entre iteraciones consecutivas:');
    disp(distances);
% Comprobación de contracción
    if all(distances(2:end) < distances(1:end-1))
        disp('La función parece ser una contracción.');
    else
        disp('La función NO parece ser una contracción.');
    end
end
