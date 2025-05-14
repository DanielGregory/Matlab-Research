# 🔬 1D Finite Difference Method – Convection-Diffusion Solver

This MATLAB project implements a **Finite Difference Method (FDM)** to solve a **singularly perturbed 1D convection-diffusion equation** with Dirichlet boundary conditions. It supports both **uniform meshes** and **Shishkin meshes** to capture boundary layers accurately and assesses numerical errors in the **maximum norm**.

---

## 📘 Problem Statement

We solve the boundary value problem:

\[
-\varepsilon u'' + u' = f(x), \quad u(0) = 0, \quad u(1) = 0
\]

where `ε` is a small positive parameter (e.g., \(10^{-6}\)) introducing boundary layer behavior.

---

## 🔧 Methodology

The code supports:

- **Uniform mesh** generation
- **(Commented)** Shishkin mesh generation for boundary-layer resolution
- Matrix assembly for diffusion and convection terms
- Computation of numerical solution using FDM
- Comparison with the **exact analytical solution**
- **Maximum norm error** and **convergence rate** evaluation

---

## 📁 File Contents

- `Research.m`: Main MATLAB script that solves the problem on varying mesh sizes, evaluates convergence, and plots results.

---

## 📊 Output

- Plots comparing **numerical vs. exact solution**
- Maximum norm error printed for each mesh size
- Convergence rates logged for increasing mesh refinements

---

## 📈 Mesh Sizes Tested

The method tests multiple uniform mesh refinements:

```matlab
p = [32, 64, 128, 256, 512, 1024, 2048, 4096, 8192];
