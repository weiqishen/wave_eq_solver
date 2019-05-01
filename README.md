# 2D Wave Equation Solver

A simple Finite Difference Solver for HPC class project.
Developers: Joshua Anton, Rishav Dutta, Hemant Kumavat, Weiqi Shen

## Building instruction

Serial Code:

1. `$g++ src2/solver2.cpp -O solver`
2. `$./solver`

CUDA Code:

1. `$nvcc src2/solver.cu -o solver`
2. `./solver`

MPI Code(Preprocessing part only):

1. `$cmake .`
2. `$ccmake .` to modify paths to libraries
3. `$make`
4. `$ ./bin/wave_eq`