# 2D Wave Equation Solver

A simple Finite Difference Solver for HPC class project.
Developers: Joshua Anton, Rishav Dutta, Hemant Kumavat, Weiqi Shen

## Code Description

This code follows methods in [Langtangen, Hans Petter. "Finite difference methods for wave motion." (2016)](http://hplgit.github.io/num-methods-for-PDEs/doc/pub/wave/pdf/wave-4print-A4.pdf) to solve the non-dimensionalized wave equation on 2D structured grid:
$$\frac{\partial ^2 u^*}{\partial t^{*2}}=\frac{\partial^2 u^*}{\partial x^{*2}}+\frac{\partial^2 u^*}{\partial y^{*2}}+f^*$$
Where $u^*=u/u_0$, $t^*=t/(L/c)$, $x^*=x/L$, $y^*=y/L$ and $f^*=f/(u_0c^2/L^2)$. Here, $u_0$ is typical speed, $L$ is length scale, $c$ is wave speed.

To utilize popular time discretization methods such as [RK34($2N^*$)](https://epubs.siam.org/doi/pdf/10.1137/07070485X), the wave equation is time splitted into two equations that are both 1st order in time.

$$ \left\{
\begin{array}{rcl}
\frac{\partial u^*}{\partial t^*}&=& v^*\\
\frac{\partial v^*}{\partial t^*}&=& \nabla^2 u^*+f^*
\end{array} \right. $$

**Initial Conditions:**
$$\begin{aligned}
u^*(x,y,0)&=&U \\
v^*(x,y,0)&=&V \\
\end{aligned} $$

The code is written in C++ language. 

**Dependencies**
[Parmetis](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview) - partition the computational mesh
[OpenMPI](https://www.open-mpi.org/) - parallelization
[Intel MKL](https://software.intel.com/en-us/mkl) - Optimized matrix and vector maths library (free for students)
[CGNS](https://cgns.github.io/) - Visualization

Below is some detailed description of the code's features

1. Code Structures
    The code cosists of the following components:
    **Mesh** - Store mesh data and repartition mesh.
    **Solver** - Initialize solution. Store solution. Calculate spatial discretization and advance in time.
    **Params** - Collect simulation parameters.
    **writer** - Write visualization file

2. Boundary conditions
    This code support 3 types of boundary conditions: 1) `Dirichlet`, 2) `Von Neumann`, 3) `Radiation`

3. Initial Condition
    For now, two types of IC is implemented. 1) `Uniform`, 2) `Gaussian pulse` 3) `Sine waves`

4. Source term
    For now, source term is hard coded to be 0.

## Building system

This project uses `CMake` as building system for convenience. To build the code, at project directory 1) `$cmake .` 2)`$ccmake .` to modify paths to libraries 3) `$make`