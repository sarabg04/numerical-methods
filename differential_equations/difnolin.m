function [X,Y,iter,incr]=difnolin(f,fy,fz,a,b,alpha,beta,N,tol,maxiter)
%DIFNOLIN  Método de diferencias finitas para problemas no lineales
%Aproxima la solución de un problema de contorno no lineal.
% RECUERDA CAMBIAR LA H (paso de discretización en x)
h = (b-a)/(N+1);
% Malla en x (incluye extremos)
X = a:h:b;
% Puntos interiores en x
x = X(2:N+1);
% Paso inicial en y (aproximación lineal entre condiciones de contorno)
k = (beta-alpha)/(N+1);
Y = alpha:k:beta;
% Valores interiores de y
y = Y(2:N+1);
% Inicialización de variables de control
incr = tol+1;   % incremento inicial mayor que tolerancia
iter = 0;       % contador de iteraciones
% Iteración principal
while incr > tol && iter < maxiter
    % Aproximación de la derivada y' en los puntos interiores
    z = (Y(3:end)-Y(1:end-2))/(2*h);
    % Evaluación de la función y sus derivadas
    fe  = feval(f,  x, y, z);
    fye = feval(fy, x, y, z);
    fze = feval(fz, x, y, z);
    % Diagonal principal del sistema
    dp = 2 + h^2*fye;
    % Subdiagonal
    ds = -1 + h/2*fze(1:end-1);
    % Superdiagonal
    di = -1 - h/2*fze(2:end);
    % Término independiente del sistema
    d = diff(Y,2) - h^2*fe;
    % Resolución del sistema tridiagonal (Crout)
    v = Crout(dp,ds,di,d);
    v = v(:);
    % Actualización de la solución
    y = y + v';
    % Reconstrucción del vector completo con condiciones de contorno
    Y = [alpha,y,beta];
    % Criterio de parada (norma del incremento)
    incr = norm(v); 
    % alternativa posible: incr = max(abs(v));
    % Incremento del contador
    iter = iter + 1;
end
end
