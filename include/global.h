#pragma once

#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <string>
#include "mpi.h"
#include "ndarray.h"

#define PI 3.14159265358979323846

struct Solution
{
    //solution values
    ndarray<ndarray<double>> sol; // local solution 2 registers
    ndarray<double> rhs;          // du/dt=rhs
};

struct Param
{
    //mesh parameters
    ndarray<double> mesh_xy0; //lower left coord
    ndarray<double> mesh_xy1; //upper right coord
    ndarray<int> mesh_nxy;    //number of points in each dir

    //simulation parameters
    double dt;
    int n_steps;
    int plot_freq;

    //non-dimensionalization
    double u_ref;
    double c_ref;
    double L_ref;
    double t_ref;

    //boundary and initial conditions
    ndarray<int> bc_type; //left right down up
    ndarray<double> bc_val;
    int ic_type, icdt_type;
    double ic_a, ic_b0, ic_b1, ic_c0, ic_c1;
    double icdt_a, icdt_b0, icdt_b1, icdt_c0, icdt_c1;

    //MPI
    int rank;
    int nproc;
};

enum BC_TYPE_ID
{
    BC_DIRICHLET,
    BC_VON_NEUMANN,
    BC_RADIATION
};

enum IC_TYPE_ID
{
    IC_UNIFORM,
    IC_SIN_WAVE,
    IC_GAUSSIAN
};