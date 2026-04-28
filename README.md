# Numerical Methods in MATLAB
This repository contains a collection of MATLAB implementations of classical and advanced numerical methods, focused on:
* Linear systems
* Nonlinear systems
* Differential equations
* Numerical integration
* Eigenvalue problems

The goal of this project is to provide a clean, well-structured, and reusable numerical methods toolkit, suitable for academic use and technical portfolios.

##  Repository Structure

```
numerical-methods-matlab/
│
├── differential_equations/
├── linear_systems/
├── nonlinear_systems/
├── numerical_integration/
├── examples/
└── README.md
```

## Linear Systems

Methods for solving systems of equations ( Ax = b ):

* **Jacobi.m**
  Classical iterative method with convergence tracking.

* **GaussSeidel.m**
  Improved iterative method with faster convergence.

* **gradiente.m**
  Gradient method for symmetric positive definite matrices.

* **gradienteconjugado.m**
  Conjugate Gradient method (efficient for large sparse systems).

* **House_Holder_QR_simetricas.m**
  QR factorization using Householder reflections.

---

## 🔹 Eigenvalue Methods

* **potenciainversa.m**
  Inverse Power Method for computing eigenvalues and eigenvectors.

---

## 🔹 Nonlinear Systems

Methods for solving ( F(x) = 0 ):

* **newton.m**
  Newton’s method for nonlinear systems.

* **jarratt.m**
  High-order multipoint iterative method.

* **OstrowskiSist.m**
  Advanced multipoint method using divided differences.

* **Stf_sistemas_vpa.m**
  Symbolic-numeric method using variable precision arithmetic (VPA).

---

## 🔹 Differential Equations

* **difnolin.m**
  Finite difference method for nonlinear boundary value problems.

---

## 🔹 Numerical Integration

* **GaussLegendre.m**
  Gaussian quadrature for double integrals.

* **Simpson.m**
  Composite Simpson’s rule for 2D integration.

* **Simpsonnorectangular.m**
  Simpson method for non-rectangular domains.

* **MonteCarlogeneral.m**
  Monte Carlo integration with geometric acceptance-rejection.

---

## 📌 Features

* Clean and modular MATLAB implementations
* Convergence tracking (residuals, increments, ACOC)
* Support for advanced methods (high-order, multipoint)
* Educational and reusable code structure

---

## ▶️ Examples

Example scripts are available in the `examples/` folder to demonstrate usage:

* `example_difnolin.m`
* `example_GaussLegendre.m`

Run them directly in MATLAB:

```matlab
run examples/example_GaussLegendre.m
```

---

## ⚠️ Notes

* Some methods require auxiliary functions such as:

  * `Crout.m`
  * `NPLegendre.m`
  * `divdiff.m`
  * `fACOC.m`

* Certain algorithms assume specific conditions:

  * Gradient methods → symmetric positive definite matrices
  * Nonlinear solvers → good initial guess required

---

## 🚀 Purpose

This repository is designed as:

* 📚 Academic support for numerical analysis courses
* 💼 Technical portfolio for computational mathematics
* 🔬 Foundation for further research or extensions

---

## 📈 Future Improvements

* Performance optimization (vectorization)
* Error analysis and benchmarking
* Visualization tools (convergence plots)
* Extended documentation

---

## 👤 Author

Developed as part of a numerical methods portfolio in MATLAB.

---

## ⭐ If you find this useful

Feel free to star ⭐ the repository or use it as a reference for your own projects.
