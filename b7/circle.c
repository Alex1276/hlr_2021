#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>


int rank, nprocs;

void printAllBuffs(int* buf)
{
    int sizeOfBuf = buf[0];
    for (int i = 0; i < nprocs; i++)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        if(rank == i)
        {   printf("(");
            for (int j = 1; j <= sizeOfBuf; j++)
            {
                printf("%d ", buf[j]);
            }
            printf(")");
        }
    }
}

void
circle(int firstElement, int* buf, int maxBufSize)
{
    MPI_Status status;
    int sendTo = rank + 1;
    int recieveFrom = rank - 1;
    int finish;

    if (rank == 0)
    {
        recieveFrom = nprocs - 1;
    }
    else if (rank == nprocs - 1)
    {
        sendTo = 0;
    }


    finish = 0;
    if (rank == nprocs - 1)
    {
        finish = (firstElement == buf[1]);
    }

    MPI_Bcast(&finish, 1, MPI_INT, nprocs - 1, MPI_COMM_WORLD);
    while(!finish)
    {
        MPI_Send(buf, maxBufSize, MPI_INT, sendTo, 1, MPI_COMM_WORLD);

        
        MPI_Recv(buf, maxBufSize, MPI_INT, recieveFrom, 1, MPI_COMM_WORLD, &status);
                
        if (rank == nprocs - 1)
        {
            finish = (firstElement == buf[1]);
        }
        MPI_Bcast(&finish, 1, MPI_INT, nprocs - 1, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
    {
        printf("\nNacher: ");
    }

    printAllBuffs(buf);
}

int*
initBuf(int maxBufSize, int bufSize)
{
    // initialize array to be the lengtg of the maximum buff size + 1
    // buf[0] will store the length of array
    int* buf = malloc(sizeof(int) * (maxBufSize + 1));

    srand(time(NULL) * (rank + 1));

    buf[0] = bufSize;
    for (int i = 1; i <= bufSize; i++)
    {
        buf[i] = rand() % 25;
    }

    return buf;
}

int main (int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int N;
    int* buf;
    int firstElement;
    
    if (argc < 2) {
        printf("Too many Arguments");
        return EXIT_FAILURE;
    }
    
    if(!sscanf(argv[1], "%i", &N))
    {
        printf("could not parse argument");
        return EXIT_FAILURE;
    }


   	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int size = N / nprocs;
    int sizeRemainder = N % nprocs;
    int bufSize;
    int maxBufSize = size + 1;

    if(rank < sizeRemainder)
    {
        bufSize = maxBufSize;
    }
    else
    {
        bufSize = size;
    }

    buf = initBuf(maxBufSize, bufSize);

    if (rank == 0)
    {
        printf("Vorher: ");
    }
    

    MPI_Barrier(MPI_COMM_WORLD);

    printAllBuffs(buf);


    firstElement = buf[1];

    MPI_Bcast(&firstElement, 1, MPI_INT, 0, MPI_COMM_WORLD);

    circle(firstElement, buf, maxBufSize);
    
    MPI_Finalize();

    return EXIT_SUCCESS;
}