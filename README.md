# 2D Wave Equation Solver

A simple Finite Difference Solver for HPC class project.
Developers: Joshua Anton, Rishav Dutta, Hemant Kumavat, Weiqi Shen

## Code Description

This code follows methods in `Langtangen, Hans Petter. "Finite difference methods for wave motion." (2016).` to solve the non-dimensionalized wave equation on 2D structured grid:
$$\frac{\partial ^2 u^*}{\partial t^{*2}}=\frac{\partial^2 u^*}{\partial x^{*2}}+\frac{\partial^2 u^*}{\partial y^{*2}}$$
Where $u^*=u_0$, $t^*=L/c$, $x^*=y^*=L$. Here, $u_0$ is typical speed, $L$ is length scale, $c$ is wave speed.

The code is written in C++ language. `Parmetis` is used to partition the grid and OpenMPI is used to parallelize the computation.

1. Code Structures
TBD
2. Boundary conditions
This code support 3 types of boundary conditions: 1) `Dirichlet`, 2) `Von Neumann`, 3) `Radiation`

3. Initial Condition
    For now, only two types of IC is implemented. 1) Gaussian pulse 2) Rectangular pulse

## Building system

This project uses `CMake` as building system for convenience. To build the code, at project directory 1) `$cmake ./` 2)`$make`