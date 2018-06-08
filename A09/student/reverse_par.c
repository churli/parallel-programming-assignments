#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h> 
#include "helper.h"
  
void reverse(char *str, int strLen)
{
	// parallelize this function and make sure to call reverse_str()
	// on each processor to reverse the substring.

	int np, rank;

    MPI_Comm_size(MPI_COMM_WORLD, &np);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int stride = strLen/np;
	int remainder = strLen%np;
    int localRemainder = 0;
    char* localBuffer = (char*) calloc(stride+1,sizeof(char)); // The +1 is to accomodate any remainder
    int* accRemainders = (int*) calloc(np,sizeof(int)); // here we put the remainders of each rank
    int* locRemainders = (int*) calloc(np,sizeof(int)); // here we put the remainders of each rank

    // If we are rank0 we must send the substrings to processes
    if (rank == 0)
    {
        // printf("Rank %d: strLen = %d\n", rank, strLen);
        int remainderToGo = remainder;
        int remainderDone = 0;
        for (int i=1; i<np; ++i)
        {
            int sendRemainder = 0;
            if (remainderToGo > 0)
            {
                sendRemainder = 1;
                locRemainders[i] = sendRemainder;
                --remainderToGo;
            }
            MPI_Send(
                str+i*stride+remainderDone,
                stride+sendRemainder,
                MPI_CHAR,
                i,
                0,
                MPI_COMM_WORLD
                );
            if (sendRemainder>0)
            {
                ++remainderDone;
            }
                accRemainders[i] = remainderDone;
        }
        if (remainderToGo > 0)
        {
            localRemainder = remainderToGo; // which should be 0, maybe to check
        }
        strncpy(localBuffer, str, stride+localRemainder);
    }
    else // If we are another rank we must receive the substring
    {
        if (rank <= remainder)
        {
            localRemainder=1;
        }
        MPI_Status status;
        MPI_Recv(
            localBuffer,
            stride+localRemainder,
            MPI_CHAR,
            0,
            0,
            MPI_COMM_WORLD,
            &status
            );
    }

    // printf("Rank %d: localBuffer = %s\n", rank, localBuffer);
    reverse_str(localBuffer, stride+localRemainder);
    // printf("Rank %d: (reverse)localBuffer = %s\n", rank, localBuffer);

    // Now we need to communicate all substrings back to rank0 and merge
    if (rank == 0)
    {
        for (int i=1; i<np; ++i)
        {
            MPI_Status status;
            int offset = strLen-((i+1)*stride+accRemainders[i]);
            int count = stride+locRemainders[i];
            // printf("Rank %d: RECV : from = %d; offset = %d; count = %d, accRemainders[%d] = %d\n", 
            // rank, i, offset, count, i, accRemainders[i]);
            MPI_Recv(
                str+offset,
                count,
                MPI_CHAR,
                i,
                MPI_ANY_TAG,
                MPI_COMM_WORLD,
                &status
                );
        }
        {
            int offset = strLen-(stride+localRemainder);
            int count = stride+localRemainder;
            // printf("Rank %d: COPY : from = %d; offset = %d; count = %d\n", rank, rank, offset, count);
            strncpy(str+offset, localBuffer, count);
        }
    }
    else // If we are another rank we must send the substring back
    {
        // printf("Rank %d: SEND : count = %d\n", rank, stride+localRemainder);
        MPI_Send(
                localBuffer,
                stride+localRemainder,
                MPI_CHAR,
                0,
                localRemainder, // we use this on the receiving end
                MPI_COMM_WORLD
                );
    }
}
