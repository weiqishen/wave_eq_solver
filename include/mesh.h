#pragma once

#include "global.h"

class mesh
{
public:

    //functions
    mesh(Param &input);
    void repartition_mesh(int nproc, int rank);
    ~mesh();

    //data
    int np_global;           //number of points global
    int np;                  //number of points local
    ndarray<double> xv;      //coordinates of local points
    ndarray<int> p2global_p; //local points id to global points id
    ndarray<int> bc_id;      //boundary condition id for each point. -1: internal; >=0 boundary; <-1 corner;
    ndarray<int> p2c;        //point to connected global point id
};
