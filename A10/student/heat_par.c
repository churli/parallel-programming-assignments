#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "heat.h"
#include "helper.h"

typedef enum Direction {TOP=0, BOTTOM=1, LEFT=2, RIGHT=3} Direction;

// We are on a rectangular grid of dimensions pX*pY
inline int getRankFromCoordinates(int pI, int pJ, int pX, int pY)
{
    if (pI < 0 || pI >= pY || pJ < 0 || pJ >= pX)
    {
        return MPI_PROC_NULL;
    }
    else
    {
        return pI*pX + pJ;
    }
}

inline void getNeighbours(int pI, int pJ, int pX, int pY, int *top, int *bot, int *left, int *right)
{
    *top = getRankFromCoordinates(pI-1, pJ, pX, pY);
    *bot = getRankFromCoordinates(pI+1, pJ, pX, pY);
    *left = getRankFromCoordinates(pI, pJ-1, pX, pY);
    *right = getRankFromCoordinates(pI, pJ+1, pX, pY);
}

double jacobi(double *h_new, double *h_old, int niters, int energy_intensity, 
            int n, int iter_energy,  const int nsources, int sources[nsources][2], 
            int rank, int size, int px, int py, MPI_Comm comm, int output)
{
    // // printf("[R%d] Allocating matrices\n", rank);
    h_old = (double *)calloc((n + 2) * (n + 2), sizeof(double)); // extended with halos of width 1
    h_new = (double *)calloc((n + 2) * (n + 2), sizeof(double)); // extended with halos of width 1

    // printf("[R%d] Allocating MPI_Request\n", rank);
    MPI_Request reqSend[4], reqRecv[4];
    // printf("[R%d] Allocating MPI_Status\n", rank);
    MPI_Status statusSend[4], statusRecv[4];
    // printf("[R%d] Declaring MPI_Datatype\n", rank);
    MPI_Datatype column;
    // printf("[R%d] Declaring MPI_Type_vector\n", rank);
    MPI_Type_vector(n, 1, n+2, MPI_DOUBLE, &column);
    // printf("[R%d] Declaring MPI_Type_commit\n", rank);
    MPI_Type_commit(&column);
    
    // Now setup MPI comms
    int pI = rank / px;
    int pJ = rank % px;
    int top, bot, left, right;
    getNeighbours(pI, pJ, px, py, &top, &bot, &left, &right);

    double *tmp;
    for (int iter = 0; iter < niters; ++iter)
    {
        for (int j = 1; j < n + 1; ++j)
        {
            for (int i = 1; i < n + 1; ++i)
            {
                h_new[map(i, j ,n+2)] = h_old[map(i, j, n+2)] / 2.0 + (h_old[map(i - 1, j, n+2)] + h_old[map(i + 1, j, n+2)] + h_old[map(i, j - 1, n+2)] + h_old[map(i, j + 1, n+2)]) / 4.0 / 2.0;
            }
        }
        if (iter < iter_energy)
        {
            for (int i = 0; i < nsources; ++i)
            {
                h_new[map(sources[i][0], sources[i][1], n+2)] += energy_intensity; // heat rate
            }
        }
        // // Debugging
        // int offset = map(0,1,n+2);
        // for (int i=0; i<n; ++i)
        // {
        //     *(h_new+offset) = 666;
        //     offset += n+2;
        // }
        // return 0.0; //debug
        // Actual MPI communication
        // printf("[R%d] Send TOP(%d)\n", rank, top);
        MPI_Isend(h_new+map(1,1,n+2), n, MPI_DOUBLE, top, 0, comm, reqSend+TOP); //Top
        // printf("[R%d] Send BOTTOM(%d)\n", rank, bot);
        MPI_Isend(h_new+map(1,n,n+2), n, MPI_DOUBLE, bot, 0, comm, reqSend+BOTTOM); //Bottom
        // printf("[R%d] Send LEFT(%d)\n", rank, left);
        MPI_Isend(h_new+map(1,1,n+2), 1, column, left, 0, comm, reqSend+LEFT); //Left
        // printf("[R%d] Send RIGHT(%d)\n", rank, right);
        MPI_Isend(h_new+map(n,1,n+2), 1, column, right, 0, comm, reqSend+RIGHT); //Right
        //
        // printf("[R%d] Recv TOP(%d)\n", rank, top);
        MPI_Irecv(h_new+map(1,0,n+2), n, MPI_DOUBLE, top, MPI_ANY_TAG, comm, reqRecv+TOP); //Top
        // printf("[R%d] Recv BOTTOM(%d)\n", rank, bot);
        MPI_Irecv(h_new+map(1,n+1,n+2), n, MPI_DOUBLE, bot, MPI_ANY_TAG, comm, reqRecv+BOTTOM); //Bottom
        // printf("[R%d] Recv LEFT(%d) -> h_new+%d\n", rank, left, map(0,1,n+2));
        MPI_Irecv(h_new+map(0,1,n+2), 1, column, left, MPI_ANY_TAG, comm, reqRecv+LEFT); //Left
        // printf("[R%d] Recv RIGHT(%d) -> h_new+%d\n", rank, right, map(n+1,1,n+2));
        MPI_Irecv(h_new+map(n+1,1,n+2), 1, column, right, MPI_ANY_TAG, comm, reqRecv+RIGHT); //Right
        //
        MPI_Waitall(4, reqSend, statusSend);
        MPI_Waitall(4, reqRecv, statusRecv);
        //
        tmp = h_new; // swap arrays
        h_new = h_old;
        h_old = tmp;
    }
    if (output) printarr(h_new, n, rank);
    return calculate_total_heat(h_new, n);
}
