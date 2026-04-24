# Numerical Methods
Colección de algoritmos numéricos implementados desde cero en MATLAB.
## Métodos incluidos
- Factorización de Cholesky
- Resolución de sistemas lineales
## Ejemplo
```matlab
A = [4 1; 1 3];
b = [1; 2];
[L, x] = cholesky_solver(A, b)


### Nonlinear Systems
- Fixed-point iteration in R²
- Contraction analysis using distance between iterates
