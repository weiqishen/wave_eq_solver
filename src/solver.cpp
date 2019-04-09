#include <algorithm>
#include "solver.h"
#include "funcs.h"

using namespace std;
solver::solver(Param *in_input, mesh *in_msh, Solution *in_solu)
{
    //store pointers
    input_ptr = in_input;
    msh_ptr = in_msh;
    solu_ptr = in_solu;

    //allocate solution/rhs array
    if (input_ptr->adv_type == EXP_EULER)
        solu_ptr->u.setup(1);
    else if (input_ptr->adv_type == RK_45)
        solu_ptr->u.setup(2);
    for (size_t i = 0; i < solu_ptr->u.get_len(); i++)
        solu_ptr->u(i).setup({(size_t)msh_ptr->np, 2});
    solu_ptr->rhs.setup({(size_t)msh_ptr->np,2});

    //setup point connectivity
    setup_connectivity();
}

void solver::setup_connectivity()
{
    //find out mpi points (in and out)
    for (size_t i = 0; i < (size_t)msh_ptr->np; i++) //loop over each local point
    {
        bool flag_mpi = false;
        for (size_t j = 0; j < 4; j++) //loop over each connection
        {
            int temp_id_global = msh_ptr->p2c({j, i});
            if (temp_id_global != -1)//not NULL(i.e at boundary)
            {
                int temp_id_local = get_index(temp_id_global, msh_ptr->p2global_p);//get local point index

                if (temp_id_local == -1) //belong to other processor, add to in_mpi_p
                {
                    flag_mpi = true;
                    in_mpi_p.push_back(temp_id_global);
                }
            }
        }
        if (flag_mpi)//if have connected point belongs to other processor
            out_mpi_p.push_back(msh_ptr->p2global_p(i)); //add to out_mpi_p
    }
    //delete repeating points and sort it in ascending order
    sort(in_mpi_p.begin(), in_mpi_p.end());
    std::vector<int>::iterator new_end = unique(in_mpi_p.begin(), in_mpi_p.end());
    in_mpi_p.resize(new_end - in_mpi_p.begin());

    //match local mpi points to other processors
    if (input_ptr->rank == 0)
        cout << "Matching MPI points... " << flush;
    match_mpi();
    if (input_ptr->rank == 0)
        cout << "Done." << endl;
    //set up pointers from solution points to buffer/neighbouring points
    if (input_ptr->rank == 0)
        cout << "Setting up pointers... " << flush;
    setup_ptrs();
    if (input_ptr->rank == 0)
        cout << "Done." << endl;
}

double *solver::get_ptr(int in_p_global, int in_field)
{
    if (in_p_global == -1)
        return NULL;
    else
    {
        int local_idx = get_index(in_p_global, msh_ptr->p2global_p);
        if (local_idx != -1) //this processor
            return solu_ptr->u(0).get_ptr({(size_t)local_idx, (size_t)in_field});
        else //in buffer
            return in_buffer.get_ptr({(size_t)in_field, (size_t)(find(in_mpi_p.begin(), in_mpi_p.end(), in_p_global) - in_mpi_p.begin())});
    }
}

void solver::send_solution()
{
    //fill out-buffer
    for (size_t i = 0; i < 2; i++)
        for (size_t j = 0; j < out_mpi_p.size(); j++)
            out_buffer({i, j}) = *ptr_out({j, i});

    //send
    size_t ct_out = 0, ct_in = 0, req_ct = 0;
    for (int p = 0; p < input_ptr->nproc; p++)
    {
        if (n_out(p) != 0)
        {
            MPI_Isend(out_buffer.get_ptr({0, ct_out}), 2 * n_out(p), MPI_DOUBLE, p, input_ptr->rank, MPI_COMM_WORLD, out_req.get_ptr(req_ct));
            MPI_Irecv(in_buffer.get_ptr({0, ct_in}), 2 * n_in(p), MPI_DOUBLE, p, p, MPI_COMM_WORLD, in_req.get_ptr(req_ct));
            ct_in += 2 * n_in(p);
            ct_out += 2 * n_out(p);
            req_ct++;
        }
    }
}

void solver::receive_solution()
{
    //wait all
    MPI_Waitall(out_req.get_len(), out_req.get_ptr(), outstatus.get_ptr());
    MPI_Waitall(in_req.get_len(), in_req.get_ptr(), outstatus.get_ptr());
}

void solver::match_mpi()
{
    int k_local_out = out_mpi_p.size(); //number of non-repeating points on this processor to send to others
    int n_req;
    ndarray<int> kproc_out(input_ptr->nproc); //number of non-repeating points on each processor to send to others
    int ct, ct1;                              //counter

    //set up n_in
    MPI_Allgather(&k_local_out, 1, MPI_INT, kproc_out.get_ptr(), 1, MPI_INT, MPI_COMM_WORLD);
    n_in.setup(input_ptr->nproc);
    vector<int> old_in_mpi_p(in_mpi_p); //copy the old one

    //first every processor broadcast out points, so that other processors can set up in points
    ct = 0;
    n_req=0;
    for (int p = 0; p < input_ptr->nproc; p++)
    {
        ndarray<int> temp_buffer(kproc_out(p));
        if (p == input_ptr->rank) //load out points
        {
            copy(out_mpi_p.begin(), out_mpi_p.end(), temp_buffer.get_ptr()); //copy to buffer
        }
        MPI_Bcast(temp_buffer.get_ptr(), kproc_out(p), MPI_INT, p, MPI_COMM_WORLD);

        if (p != input_ptr->rank) //other processors are broadcasting
        {
            vector<int> intersect;
            set_intersection(temp_buffer.get_ptr(), temp_buffer.get_ptr() + kproc_out(p),
                             old_in_mpi_p.begin(), old_in_mpi_p.end(), back_inserter(intersect));
            n_in(p) = intersect.size();
            if (n_in(p) != 0) //have points to receive from p
            {
                n_req++;
                copy(intersect.begin(), intersect.end(), in_mpi_p.data() + ct); //the size of vector will not change due to non-repeating nature
                ct += n_in(p);
            }
        }
        else
        {
            n_in(p) = 0; //nothing to receive from self
        }
    }

    //calculate n_out
    n_out.setup(input_ptr->nproc);
    MPI_Alltoall(n_in.get_ptr(), 1, MPI_INT, n_out.get_ptr(), 1, MPI_INT, MPI_COMM_WORLD); //notice that n_out(rank)=0
    size_t temp_size = 0;                                                                  //total number of sending points(allow repeating)
    for (size_t i = 0; i < n_out.get_len(); i++)
        temp_size += n_out(i);
    out_mpi_p.resize(temp_size); //change size of outmpi_p due to repeating

    //allocate memory for requests
    in_req.setup(n_req);
    out_req.setup(n_req);
    instatus.setup(n_req);
    outstatus.setup(n_req);

    //fill out_mpi_p
    ct = 0, ct1 = 0;
    size_t req_ct = 0;
    for (int p = 0; p < input_ptr->nproc; p++) //loop over each proc
    {
        if (n_in(p) != 0) //if need to send in_mpi to this processor
        {
            MPI_Isend(in_mpi_p.data() + ct, n_in(p), MPI_INT, p, input_ptr->rank, MPI_COMM_WORLD, in_req.get_ptr(req_ct));
            MPI_Irecv(out_mpi_p.data() + ct1, n_out(p), MPI_INT, p, p, MPI_COMM_WORLD, out_req.get_ptr(req_ct));
            ct += n_in(p);
            ct1 += n_out(p);
            req_ct++;
        }
    }
    MPI_Waitall(n_req, in_req.get_ptr(), instatus.get_ptr());
    MPI_Waitall(n_req, out_req.get_ptr(), outstatus.get_ptr());

    //set up in and out buffer
    in_buffer.setup({2, in_mpi_p.size()});
    out_buffer.setup({2, out_mpi_p.size()});
}

void solver::setup_ptrs()
{
    //set up out sendinng pointers
    ptr_out.setup({out_mpi_p.size(), 2});
    for (size_t i = 0; i < 2; i++)
        for (size_t j = 0; j < out_mpi_p.size(); j++)
            ptr_out({j, i}) = get_ptr(out_mpi_p[j], i);

    //set up pointers to neighbour of local points
    ptr_ngb.setup({4, (size_t)msh_ptr->np, 2});
    for (size_t i = 0; i < 2; i++)//field
        for (size_t j = 0; j < (size_t)msh_ptr->np; j++)//point
            for (size_t k = 0; k < 4; k++)//neighbours
                ptr_ngb({k, j, i}) = get_ptr(msh_ptr->p2c({k, j}), i);
}

void solver::set_ic()
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
