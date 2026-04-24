function I = MonteCarlogeneral(f, n, a, b, c, d)
% MONTECARLOGENERAL
% Estimación de integrales mediante método de Monte Carlo
% en regiones bidimensionales con restricciones geométricas
% INPUT:
%   f -> función f(x,y)
%   n -> número de puntos aleatorios
%   [a,b] -> intervalo en x
%   [c,d] -> intervalo en y
% OUTPUT:
%   I -> aproximación de la integral
k = 1;
% Puntos aceptados
x = [];
y = [];
% Generación de puntos aleatorios
while k <= n
    % Punto aleatorio en el dominio [a,b] x [c,d]
    u = rand;
    u = a + (b - a) * u;
    v = rand;
    v = c + (d - c) * v;
    % Condición geométrica (región de integración)
    % IMPORTANTE: definir borde1 y borde2 según el problema
    w = borde1 - v;
    t = borde2 - v;
    % Rechazo / aceptación del punto
    if w <= 0 && t >= 0
        x(k) = u;
        y(k) = v;
        k = k + 1;
    end
end
% Visualización de puntos aceptados
plot(x, y, '.')
grid on
% Área de la región
% (depende del problema: debe definirse correctamente)
[Area, ~] = simpson(@(x) resta_de_bordes(x), -1, 1, 20);
% Evaluación de la función
fi = feval(f, x, y);
% Estimación Monte Carlo
I = (1/n) * Area * sum(fi);
end
