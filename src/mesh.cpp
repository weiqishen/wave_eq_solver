
#include "mesh.h"
#include "metis.h"
#include "parmetis.h"

using namespace std;

mesh::mesh(Param &input)
{
    int start;
    double dx, dy;
    int x_i, y_i;

    //generating partial mesh in naive order
    if (input.rank == 0)
        cout << "Generating mesh connectivity... " << flush;
    np_global = input.mesh_nxy(0) * input.mesh_nxy(1);
    np = np_global / input.nproc;
    start = np * input.rank;
    if (input.rank == input.nproc - 1)
        np += np_global - input.nproc * np;

    //setup points
    p2global_p.setup(np);
    p2c.setup({4, (size_t)np});
    p2c = -1; //initialize to -1, which means NULL
    bc_id.setup(np);

    for (size_t i = 0; i < (size_t)np; i++)
    {
        p2global_p(i) = (int)i + start;
        x_i = p2global_p(i) % input.mesh_nxy(0);
        y_i = p2global_p(i) / input.mesh_nxy(0);

        //corners
        if (x_i == 0 && y_i == 0) //bottom left
        {
            bc_id(i) = -2;
            p2c({0, i}) = p2global_p(i) + 1;
            p2c({1, i}) = p2global_p(i) + input.mesh_nxy(0);
        }
        else if (x_i == input.mesh_nxy(0) - 1 && y_i == 0) //bottom right
        {
            bc_id(i) = -3;
            p2c({0, i}) = p2global_p(i) - 1;
            p2c({1, i}) = p2global_p(i) + input.mesh_nxy(0);
        }
        else if (x_i == 0 && y_i == input.mesh_nxy(1) - 1) //top left
        {
            bc_id(i) = -4;
            p2c({0, i}) = p2global_p(i) + 1;
            p2c({1, i}) = p2global_p(i) - input.mesh_nxy(0);
        }
        else if (x_i == input.mesh_nxy(0) - 1 && y_i == input.mesh_nxy(1) - 1) //top right
        {
            bc_id(i) = -5;
            p2c({0, i}) = p2global_p(i) - 1;
            p2c({1, i}) = p2global_p(i) - input.mesh_nxy(0);
        }
        else
        {
            //boundaries
            if (x_i == 0) //left middle
            {
                bc_id(i) = 0;
                p2c({0, i}) = p2global_p(i) + 1;
                p2c({1, i}) = p2global_p(i) - input.mesh_nxy(0);
                p2c({2, i}) = p2global_p(i) + input.mesh_nxy(0);
            }
            else if (x_i == input.mesh_nxy(0) - 1) //right middle
            {
                bc_id(i) = 1;
                p2c({0, i}) = p2global_p(i) - 1;
                p2c({1, i}) = p2global_p(i) - input.mesh_nxy(0);
                p2c({2, i}) = p2global_p(i) + input.mesh_nxy(0);
            }
            else if (y_i == 0) //bottom middle
            {
                bc_id(i) = 2;
                p2c({0, i}) = p2global_p(i) + input.mesh_nxy(0);
                p2c({1, i}) = p2global_p(i) - 1;
                p2c({2, i}) = p2global_p(i) + 1;
            }
            else if (y_i == input.mesh_nxy(1) - 1) //top middle
            {
                bc_id(i) = 3;
                p2c({0, i}) = p2global_p(i) - input.mesh_nxy(0);
                p2c({1, i}) = p2global_p(i) - 1;
                p2c({2, i}) = p2global_p(i) + 1;
            }
            else //internal
            {
                bc_id(i) = -1;
                p2c({0, i}) = p2global_p(i) - 1;
                p2c({1, i}) = p2global_p(i) + 1;
                p2c({2, i}) = p2global_p(i) - input.mesh_nxy(0);
                p2c({3, i}) = p2global_p(i) + input.mesh_nxy(0);
            }
        }
    }
    if (input.rank == 0)
        cout << "Done." << endl;

    repartition_mesh(input.nproc, input.rank);

    //calculate coordinates
    if (input.rank == 0)
        cout << "Calculating mesh point coordinates... " << flush;
    xv.setup({2, (size_t)np});
    dx = (double)(input.mesh_xy1(0) - input.mesh_xy0(0)) / (input.mesh_nxy(0) - 1.);
    dy = (double)(input.mesh_xy1(1) - input.mesh_xy0(1)) / (input.mesh_nxy(1) - 1.);
    for (size_t i = 0; i < (size_t)np; i++)
    {
        x_i = p2global_p(i) % input.mesh_nxy(0);
        y_i = p2global_p(i) / input.mesh_nxy(0);
        xv({0, i}) = input.mesh_xy0(0) + dx * x_i;
        xv({1, i}) = input.mesh_xy0(1) + dy * y_i;
    }
    if (input.rank == 0)
        cout << "Done." << endl;
}

void mesh::repartition_mesh(int nproc, int rank)
{
    idx_t wgtflag = 0; //no weighting
    idx_t numflag = 0; //C style
    idx_t ncon = 1;    //only on constrant
    idx_t edgecut;
    idx_t options[3];
    size_t ct;           //counter
    real_t ubvec = 1.05; //imbalance tolerance
    ndarray<idx_t> xadj, adjncy, vtxdist;
    ndarray<idx_t> kproc(nproc);   //number of points per processor
    ndarray<real_t> tpwgts(nproc); //fraction of weight per processor
    ndarray<idx_t> part(np);
    MPI_Comm comm;

    MPI_Allgather(&np, 1, MPI_INT, kproc.get_ptr(), 1, MPI_INT, MPI_COMM_WORLD);
    vtxdist.setup(nproc + 1);
    vtxdist(0) = 0;
    for (size_t i = 0; i < (size_t)nproc; i++)
        vtxdist(i + 1) = vtxdist(i) + kproc(i);

    xadj.setup(np + 1);
    xadj(0) = 0;
    for (size_t i = 0; i < (size_t)np; i++)
    {
        if (bc_id(i) < -1) //corner
            xadj(i + 1) = xadj(i) + 2;
        else if (bc_id(i) == -1) //internal
            xadj(i + 1) = xadj(i) + 4;
        else //boundary
            xadj(i + 1) = xadj(i) + 3;
    }

    adjncy.setup(xadj(np));
    ct = 0;
    for (size_t i = 0; i < (size_t)np; i++)                          //for each node
        for (size_t j = 0; j < (size_t)(xadj(i + 1) - xadj(i)); j++) //for each connection
            adjncy(ct++) = p2c({j, i});

    MPI_Comm_dup(MPI_COMM_WORLD, &comm);
    options[0] = 1;
    options[1] = 7;
    options[2] = 0;
    tpwgts = 1. / nproc;

    if (rank == 0)
        cout << "Start partitioning mesh... " << endl;
    if (ParMETIS_V3_PartKway(vtxdist.get_ptr(), xadj.get_ptr(), adjncy.get_ptr(),
                             NULL, NULL, &wgtflag,
                             &numflag, &ncon, &nproc,
                             tpwgts.get_ptr(), &ubvec, options,
                             &edgecut, part.get_ptr(), &comm) != METIS_OK)
        Fatal_Error("Failed to partition mesh");

    // Now create sending buffer
    ndarray<int> outK(nproc); //on this processor, number of points to each processor
    outK = 0;
    for (size_t i = 0; i < (size_t)np; i++)
        ++outK(part(i));

    ndarray<int> inK(nproc); //number of points to receive from each processor

    MPI_Alltoall(outK.get_ptr(), 1, MPI_INT, inK.get_ptr(), 1, MPI_INT, MPI_COMM_WORLD);

    ndarray<int> newkprocs(nproc); //total number of points on each processor
    MPI_Allreduce(outK.get_ptr(), newkprocs.get_ptr(), nproc, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    int totalinK = newkprocs(rank); //total number of points belong to this processor

    //number of reveive/send requests
    int num_inrequests = 0;
    int num_outrequests = 0;
    for (size_t i = 0; i < (size_t)nproc; i++)
    {
        if (inK(i) != 0)
            num_inrequests++;
        if (outK(i) != 0)
            num_outrequests++;
    }

    //declare sending buffer
    ndarray<ndarray<int>> outp2c(num_outrequests);        //p2c send to each processor
    ndarray<ndarray<int>> outbc_id(num_outrequests);      //bc_id send to each processor
    ndarray<ndarray<int>> outp2global_p(num_outrequests); //p2global_p send to each processor
    // declare recieving buffer
    ndarray<int> new_p2c({4, (size_t)totalinK});
    ndarray<int> new_bc_id(totalinK);
    ndarray<int> new_p2global_p(totalinK);

    //declare mpi requests and status
    ndarray<MPI_Request> inrequests_p2c(num_inrequests);
    ndarray<MPI_Request> inrequests_bc_id(num_inrequests);
    ndarray<MPI_Request> inrequests_p2global_p(num_inrequests);

    ndarray<MPI_Request> outrequests_p2c(num_outrequests);
    ndarray<MPI_Request> outrequests_bc_id(num_outrequests);
    ndarray<MPI_Request> outrequests_p2global_p(num_outrequests);

    ndarray<MPI_Status> instatus(num_inrequests);
    ndarray<MPI_Status> outstatus(num_outrequests);

    // Make exchange for arrays

    ct = 0;
    size_t inrequest_counter = 0;
    for (int p = 0; p < nproc; p++) //for each processor, ensure ascending order of p2global_p
    {
        if (inK(p) != 0) //if have data to receive
        {
            MPI_Irecv(new_p2c.get_ptr({0, ct}), 4 * inK(p), MPI_INT, p, p, MPI_COMM_WORLD, inrequests_p2c.get_ptr(inrequest_counter));
            MPI_Irecv(new_bc_id.get_ptr(ct), inK(p), MPI_INT, p, nproc + p, MPI_COMM_WORLD, inrequests_bc_id.get_ptr(inrequest_counter));
            MPI_Irecv(new_p2global_p.get_ptr(ct), inK(p), MPI_INT, p, 2 * nproc + p, MPI_COMM_WORLD, inrequests_p2global_p.get_ptr(inrequest_counter));
            ct += inK(p);
            inrequest_counter++;
        }
    }

    size_t outrequest_counter = 0;
    for (int p = 0; p < nproc; p++)
    {
        if (outK(p) != 0) //if this processor have points to send to processor p
        {
            ct = 0;

            outp2c(outrequest_counter).setup({4, (size_t)outK(p)});
            outbc_id(outrequest_counter).setup(outK(p));
            outp2global_p(outrequest_counter).setup(outK(p));

            for (size_t i = 0; i < (size_t)np; i++) //loop over all local elements, ensure ascending order of p2global_p
            {
                if (part(i) == p) //if this element send to processor p
                {
                    for (size_t v = 0; v < 4; v++)
                        outp2c(outrequest_counter)({v, ct}) = p2c({v, i});

                    outbc_id(outrequest_counter)(ct) = bc_id(i);
                    outp2global_p(outrequest_counter)(ct) = p2global_p(i);
                    ct++;
                }
            }
            MPI_Isend(outp2c(outrequest_counter).get_ptr(), 4 * outK(p), MPI_INT, p, rank, MPI_COMM_WORLD, outrequests_p2c.get_ptr(outrequest_counter));
            MPI_Isend(outbc_id(outrequest_counter).get_ptr(), outK(p), MPI_INT, p, nproc + rank, MPI_COMM_WORLD, outrequests_bc_id.get_ptr(outrequest_counter));
            MPI_Isend(outp2global_p(outrequest_counter).get_ptr(), outK(p), MPI_INT, p, 2 * nproc + rank, MPI_COMM_WORLD, outrequests_p2global_p.get_ptr(outrequest_counter));
            outrequest_counter++;
        }
    }

    MPI_Waitall(num_inrequests, inrequests_p2c.get_ptr(), instatus.get_ptr());
    MPI_Waitall(num_inrequests, inrequests_bc_id.get_ptr(), instatus.get_ptr());
    MPI_Waitall(num_inrequests, inrequests_p2global_p.get_ptr(), instatus.get_ptr());

    MPI_Waitall(num_outrequests, outrequests_p2c.get_ptr(), outstatus.get_ptr());
    MPI_Waitall(num_outrequests, outrequests_bc_id.get_ptr(), outstatus.get_ptr());
    MPI_Waitall(num_outrequests, outrequests_p2global_p.get_ptr(), outstatus.get_ptr());

    //copy back
    p2c = new_p2c;
    bc_id = new_bc_id;
    p2global_p = new_p2global_p; //increasing order

    np = totalinK; //update new local cell number

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        cout << "Done partitioning mesh" << endl;
}

mesh::~mesh()
{
}
