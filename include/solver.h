#pragma once

#include "global.h"
#include "mesh.h"
class solver
{
  public:
    solver(Param *in_input,mesh *in_msh,Solution *in_solu);
    void set_ic();//set initial condition
    void calc_spatial();//calculate spatial discretization (rhs)
    void time_advance();//time integration
    ~solver();

  private:
  void setup_connectivity();
  void match_mpi();
  void setup_ptrs();//set up pointers to solution array or buffers
  double* get_ptr_u(int in_p_global);

  void send_solution();
  void receive_solution();

  //ptrs
  Param *input_ptr;
  mesh *msh_ptr;
  Solution *solu_ptr;
  ndarray<double *> ptr_ngb;

  //MPI
  ndarray<double *> ptr_out;//pointers to the solution that should be send to other processors
  ndarray<double> out_buffer, in_buffer; //sending/recievinng buffer
  vector<int> in_mpi_p, out_mpi_p;//list of global point index that send to/receive from each processor
  ndarray<int> n_in, n_out;//number of points need to send to/receive from each processor
  ndarray<MPI_Request> in_req, out_req;
  ndarray<MPI_Status> instatus, outstatus;
};
