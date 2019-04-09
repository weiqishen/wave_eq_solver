
#include "reader.h"
#include "param_reader.h"

using namespace std;

void read_params(char *input_fnameC, Param &input)
{
    string input_fnameS(input_fnameC);
    param_reader pr(input_fnameS);

    if (input.rank == 0)
        cout << "Reading input parameters... " << flush;

    pr.openFile();
    //mesh parameters
    pr.getVectorValue("mesh_xy0", input.mesh_xy0);
    pr.getVectorValue("mesh_xy1", input.mesh_xy1);
    pr.getVectorValue("mesh_nxy", input.mesh_nxy);
    if (input.mesh_nxy.get_len() != 2 || input.mesh_xy0.get_len() != 2 || input.mesh_xy1.get_len() != 2)
        Fatal_Error("Unexpected number of mesh parameters");
    //simulation parameters
    pr.getScalarValue("dt", input.dt);
    pr.getScalarValue("adv_type", input.adv_type);
    pr.getScalarValue("n_steps", input.n_steps);
    pr.getScalarValue("plot_freq", input.plot_freq);

    //non-dimensionalization
    pr.getScalarValue("u_ref", input.u_ref);
    pr.getScalarValue("c_ref", input.c_ref);
    pr.getScalarValue("L_ref", input.L_ref);

    //boundary and initial conditions
    pr.getVectorValue("bc_type", input.bc_type); //left right down up
    pr.getVectorValue("bc_val", input.bc_val);   //left right down up
    if (input.bc_type.get_len() != 4 || input.bc_val.get_len() != 4)
        Fatal_Error("Number of boundary conditions not equal to 4");
    pr.getScalarValue("ic_type", input.ic_type);
    pr.getScalarValue("icdt_type", input.icdt_type);

    pr.getScalarValue("ic_a", input.ic_a);
    pr.getScalarValue("icdt_a", input.icdt_a);

    if (input.ic_type != IC_UNIFORM)
    {
        pr.getScalarValue("ic_b0", input.ic_b0);
        pr.getScalarValue("ic_b1", input.ic_b1);
        pr.getScalarValue("ic_c0", input.ic_c0);
        pr.getScalarValue("ic_c1", input.ic_c1);
    }
    if (input.icdt_type != IC_UNIFORM)
    {
        pr.getScalarValue("icdt_b0", input.icdt_b0);
        pr.getScalarValue("icdt_b1", input.icdt_b1);
        pr.getScalarValue("icdt_c0", input.icdt_c0);
        pr.getScalarValue("icdt_c1", input.icdt_c1);
    }

    pr.closeFile();
    if (input.rank == 0)
        cout << "Done." << endl;



    //non-dimensionalization
    if (input.rank == 0)
        cout << "Non-dimensionalizing... " << flush;
    //mesh parameters
    for (size_t i = 0; i < 2; i++)
    {
        input.mesh_xy0(i) /= input.L_ref; //lower left coord
        input.mesh_xy1(i) /= input.L_ref; //upper right coord
    }

    //simulation parameters
    input.t_ref = input.L_ref / input.c_ref;
    input.dt /= input.t_ref;

    //boundary conditions
    for (size_t i = 0; i < 4; i++)
    {
        if (input.bc_type(i) == BC_DIRICHLET)
            input.bc_val(i) /= input.u_ref;
        else if (input.bc_type(i) == BC_VON_NEUMANN)
            input.bc_val(i) /= (input.u_ref / input.L_ref);
        else if (input.bc_type(i) == BC_RADIATION)
            input.bc_val(i) /= (input.u_ref / input.t_ref);
        else
            Fatal_Error("Unsupported boundary condition")
    }

    //initial condition
    if (input.ic_type != IC_UNIFORM && input.ic_type != IC_SIN_WAVE && input.ic_type != IC_GAUSSIAN)
        Fatal_Error("Unsupported initial condition u_0");
    if (input.icdt_type != IC_UNIFORM && input.icdt_type != IC_SIN_WAVE && input.icdt_type != IC_GAUSSIAN)
        Fatal_Error("Unsupported initial condition u_t0");

    input.ic_a /= input.u_ref;

    if (input.ic_type == IC_GAUSSIAN)
    {
        //Gaussian: u=a*e^(-(x-b0)^2/(2*c0^2)-(y-b1)^2/(2*c1^2))
        input.ic_b0 /= input.L_ref;
        input.ic_b1 /= input.L_ref;
        input.ic_c0 /= input.L_ref;
        input.ic_c1 /= input.L_ref;
    }
    else if (input.ic_type == IC_SIN_WAVE)
    {
        //sin wave: u=a*sin(b0*(x-c0))sin(b1*(y-c1))
        input.ic_b0 *= input.L_ref;
        input.ic_b1 *= input.L_ref;
        input.ic_c0 /= input.L_ref;
        input.ic_c1 /= input.L_ref;
    }

    if (input.icdt_type == IC_GAUSSIAN)
    {
        //Gaussian: u_t=a*e^(-(x-b0)^2/(2*c0^2)-(y-b1)^2/(2*c1^2))
        input.icdt_b0 /= input.L_ref;
        input.icdt_b1 /= input.L_ref;
        input.icdt_c0 /= input.L_ref;
        input.icdt_c1 /= input.L_ref;
    }
    else if (input.icdt_type == IC_SIN_WAVE)
    {
        //sin wave: u_t=a*sin(b0*(x-c0))sin(b1*(y-c1))
        input.icdt_b0 *= input.L_ref;
        input.icdt_b1 *= input.L_ref;
        input.icdt_c0 /= input.L_ref;
        input.icdt_c1 /= input.L_ref;
    }
    
    if (input.rank == 0)
        cout << "Done." << endl;
}
