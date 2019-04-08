#include "solver.h"

using namespace std;
solver::solver(Param &in_input, mesh &in_msh)
{
    input_ptr = &in_input;
    msh_ptr = &in_msh;
}

void solver::initialize_solution()
{

}

void solver::calc_spatial()
{
}
void solver::time_advance()
{
}
solver::~solver()
{
}
