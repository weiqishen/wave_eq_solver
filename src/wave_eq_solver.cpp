#include "global.h"
#include "solver.h"
#include "reader.h"
#include "mesh.h"

using namespace std;

int main(int argc, char *argv[])
{
    int i_steps = 0;
    Solution wave_solution;
    Param input;

    //initialize MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &input.rank);
    MPI_Comm_size(MPI_COMM_WORLD, &input.nproc);
    //read args
    if (argc < 2)
    {
        if (input.rank == 0)
            cout << "No input file specified!" << endl;
        MPI_Finalize();
        return 0;
    }
    else
    {
        if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help"))
        {
            if (input.rank == 0)
            {
                cout << "This code is a simple 2D finite difference wave equation solver" << endl;
                cout << "https://github.com/weiqishen/wave_eq_solver" << endl;
            }
            MPI_Finalize();
            return 0;
        }
    }

    //read input params
    read_params(argv[1], input);
    //generate mesh
    mesh msh(input);
    //initialize solver and solution
    solver solver_wave(input, msh);
    solver_wave.initialize_solution();

    //main loop
    while (i_steps < input.n_steps)
    {
        //solve

        i_steps++;
        //write result
        if (i_steps % input.plot_freq)
        {
        }
    }

    MPI_Finalize();

    return 0;
}