#pragma once

#include "global.h"
#include "mesh.h"
class solver
{
  public:
    solver(Param &in_input,mesh &in_msh);

    void initialize_solution();
    void calc_spatial();
    void time_advance();
    ~solver();

  private:
    Param *input_ptr;
    mesh *msh_ptr;
    ndarray<double> out_buffer,in_buffer;
    //solution pointers
    ndarray<double *> p2nb,out2p; //pointers to neighbouring points (left right down up), pointers to out buffer

};
