#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "heat.h"
#include "helper.h"

typedef enum Direction {TOP=0, BOTTOM=1, LEFT=2, RIGHT=3} Direction;

static inline int mapp(const int i,const int j,const int n){return j*n+i;}

// We are on a rectangular grid of dimensions pX*pY
static inline int getRankFromCoordinates(int pI, int pJ, int pX, int pY)
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

static inline void getNeighbours(int pI, int pJ, int pX, int pY, int *top, int *bot, int *left, int *right)
{
    *top = getRankFromCoordinates(pI-1, pJ, pX, pY);
    *bot = getRankFromCoordinates(pI+1, pJ, pX, pY);
    *left = getRankFromCoordinates(pI, pJ-1, pX, pY);
    *right = getRankFromCoordinates(pI, pJ+1, pX, pY);
}

static double local_calculate_total_heat_omp(double *h_new, int numRows, int numColumns, int num_threads)
{
    double heat = 0.0; // total heat in system
    #pragma omp parallel for num_threads(num_threads) reduction(+:heat)
    for (int i = 1; i < numColumns + 1; ++i)
    {
      for (int j = 1; j < numRows + 1; ++j)
      {
        heat += h_new[mapp(i, j, numColumns+2)];
      }
    }
    return heat;
}

double jacobi(double *h_new, double *h_old, int niters, int energy_intensity, 
  int n, int iter_energy, const int nsources, int sources[nsources][2], 
  int rank, int size, int px, int py, int num_threads, MPI_Comm comm, int output)
{
    int numRows = n/py;
    int numColumns = n/px;

    h_old = (double *)calloc((numRows + 2) * (numColumns + 2), sizeof(double)); // extended with halos of width 1
    h_new = (double *)calloc((numRows + 2) * (numColumns + 2), sizeof(double)); // extended with halos of width 1

    MPI_Request reqSend[4], reqRecv[4];
    MPI_Status statusSend[4], statusRecv[4];
    MPI_Datatype column;
    MPI_Type_vector(numRows, 1, numColumns+2, MPI_DOUBLE, &column);
    MPI_Type_commit(&column);
    
    // Now setup MPI comms
    int pI = rank / px;
    int pJ = rank % px;
    int top, bot, left, right;
    getNeighbours(pI, pJ, px, py, &top, &bot, &left, &right);

    int localColRangeStart = numColumns*pJ;
    int localColRangeEnd = numColumns*(pJ+1);
    int localRowRangeStart = numRows*pI;
    int localRowRangeEnd = numRows*(pI+1);

    double *tmp;
    int commuStride = 10;
    for (int iter = 0; iter < niters; iter+=commuStride)
    {
        // Launching the non-blocking receive operations
        MPI_Irecv(h_new+mapp(1,0,numColumns+2), numColumns, MPI_DOUBLE, top, MPI_ANY_TAG, comm, reqRecv+TOP); //Top
        MPI_Irecv(h_new+mapp(1,numRows+1,numColumns+2), numColumns, MPI_DOUBLE, bot, MPI_ANY_TAG, comm, reqRecv+BOTTOM); //Bottom
        MPI_Irecv(h_new+mapp(0,1,numColumns+2), 1, column, left, MPI_ANY_TAG, comm, reqRecv+LEFT); //Left
        MPI_Irecv(h_new+mapp(numColumns+1,1,numColumns+2), 1, column, right, MPI_ANY_TAG, comm, reqRecv+RIGHT); //Right
        //
        int locIterTarget = iter + commuStride - 1;
        // for (int locIter = 0; locIter < locIterTarget ; ++locIter)
        // {
        //     for (int j = 1; j < numRows + 1; ++j)
        //     {
        //         for (int i = 1; i < numColumns + 1; ++i)
        //         {
        //             h_new[mapp(i, j ,numColumns+2)] = h_old[mapp(i, j, numColumns+2)] / 2.0 + (h_old[mapp(i - 1, j, numColumns+2)] + h_old[mapp(i + 1, j, numColumns+2)] + h_old[mapp(i, j - 1, numColumns+2)] + h_old[mapp(i, j + 1, numColumns+2)]) / 4.0 / 2.0;
        //         }
        //     }
        //     if (iter+locIter < iter_energy)
        //     {
        //         for (int i = 0; i < nsources; ++i)
        //         {
        //             int row = sources[i][1];
        //             int col = sources[i][0];
        //             if (localRowRangeStart <= row && row < localRowRangeEnd 
        //                 && localColRangeStart <= col && col < localColRangeEnd)
        //             {
        //                 row = row % numRows;
        //                 col = col % numColumns;
        //                 h_new[mapp(col, row, numColumns+2)] += energy_intensity; // heat rate
        //             }
        //         }
        //     }
        //     tmp = h_new; // swap arrays
        //     h_new = h_old;
        //     h_old = tmp;
        // }

        // #pragma omp parallel for num_threads(num_threads)
        for (int i = 1; i < numColumns + 1; ++i)
        {
            int j = 1;
            h_new[mapp(i, j ,numColumns+2)] = h_old[mapp(i, j, numColumns+2)] / 2.0 + (h_old[mapp(i - 1, j, numColumns+2)] + h_old[mapp(i + 1, j, numColumns+2)] + h_old[mapp(i, j - 1, numColumns+2)] + h_old[mapp(i, j + 1, numColumns+2)]) / 4.0 / 2.0;
        }
        MPI_Isend(h_new+mapp(1,1,numColumns+2), numColumns, MPI_DOUBLE, top, 0, comm, reqSend+TOP); //Top
        // #pragma omp parallel for num_threads(num_threads)
        for (int j = 2; j < numRows; ++j)
        {
            for (int i = 1; i < numColumns + 1; ++i)
            {
                h_new[mapp(i, j ,numColumns+2)] = h_old[mapp(i, j, numColumns+2)] / 2.0 + (h_old[mapp(i - 1, j, numColumns+2)] + h_old[mapp(i + 1, j, numColumns+2)] + h_old[mapp(i, j - 1, numColumns+2)] + h_old[mapp(i, j + 1, numColumns+2)]) / 4.0 / 2.0;
            }
        }
        MPI_Isend(h_new+mapp(1,1,numColumns+2), 1, column, left, 0, comm, reqSend+LEFT); //Left
        MPI_Isend(h_new+mapp(numColumns,1,numColumns+2), 1, column, right, 0, comm, reqSend+RIGHT); //Right
        // #pragma omp parallel for num_threads(num_threads)
        for (int i = 1; i < numColumns + 1; ++i)
        {
            int j = numRows;
            h_new[mapp(i, j ,numColumns+2)] = h_old[mapp(i, j, numColumns+2)] / 2.0 + (h_old[mapp(i - 1, j, numColumns+2)] + h_old[mapp(i + 1, j, numColumns+2)] + h_old[mapp(i, j - 1, numColumns+2)] + h_old[mapp(i, j + 1, numColumns+2)]) / 4.0 / 2.0;
        }
        MPI_Isend(h_new+mapp(1,numRows,numColumns+2), numColumns, MPI_DOUBLE, bot, 0, comm, reqSend+BOTTOM); //Bottom
        if (iter+locIterTarget < iter_energy)
        {
            for (int i = 0; i < nsources; ++i)
            {
                int row = sources[i][1];
                int col = sources[i][0];
                if (localRowRangeStart <= row && row < localRowRangeEnd 
                    && localColRangeStart <= col && col < localColRangeEnd)
                {
                    row = row % numRows;
                    col = col % numColumns;
                    h_new[mapp(col, row, numColumns+2)] += energy_intensity; // heat rate
                }
            }
        }
        //
        MPI_Waitall(4, reqSend, statusSend);
        MPI_Waitall(4, reqRecv, statusRecv);
        //
        tmp = h_new; // swap arrays
        h_new = h_old;
        h_old = tmp;
    }
    if (output) printarr(h_new, numRows, rank);

    return local_calculate_total_heat_omp(h_new, numRows, numColumns, num_threads);
}
