#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <mpi.h>


int rank, nprocs;

// function that prints all buffs from all processes
void printAllBuffs(int* buf, int sizeOfBuf)
{
    // for each process
    for (int i = 0; i < nprocs; i++)
    {
      // wait, so its in the right order
        MPI_Barrier(MPI_COMM_WORLD);
        // if process is rank that should print right now
        if(rank == i)
        {   // print all values in buffer
          printf("(");

            for (int j = 0; j < sizeOfBuf; j++)
            {
                printf("%d ", buf[j]);
            }
            printf(")");
        }
    }
}


// circle methode rotates array
void
circle(int firstElement, int* buf, int bufSize)
{
    // status to use in probe
    MPI_Status status;
    // the rank we want to send the buffer to
    int sendTo = (rank + 1) % nprocs;
    // the rank we want to receive a buffer from
    int receiveFrom;
    if (rank - 1 == -1)
    {
      receiveFrom = nprocs - 1;
    }
    else
    {
      receiveFrom = rank -1;
    }
    // 1 = done, 0 = not done yet
    int finish = 0;

    // while not finished
    while(!finish)
    {
      // send buf to process with rank sendTo
        MPI_Send(buf, bufSize, MPI_INT, sendTo, 1, MPI_COMM_WORLD);

        // get nect buffer size
        MPI_Probe(receiveFrom, 1, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &bufSize);

        buf = realloc(buf, sizeof(int) * bufSize);
        
        // receive new buffer in buf
        MPI_Recv(buf, bufSize, MPI_INT, receiveFrom, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        
        // if this is process with last rank
        if (rank == nprocs - 1)
        {   
          // done if first firstElement is first element of our buffer
            finish = (firstElement == buf[0]);
        }
        // send finished to all processes from process with last rank
        MPI_Bcast(&finish, 1, MPI_INT, nprocs - 1, MPI_COMM_WORLD);
    }

    // Done 
    // rank 0 prints Nachher
    if (rank == 0)
    {
        printf("\nNacher: ");
    }
    // all processes print their buf in right order
    printAllBuffs(buf, bufSize);
}

// allocates buffer on heap with desired size and fills random int from 0 to 24
int*
initBuf(int bufSize)
{
    //allocates memory of desired size
    int* buf = malloc(sizeof(int) * bufSize);
    // set random seed
    srand(time(NULL) * (rank + 1));
    // init values into buffer
    for (int i = 0; i < bufSize; i++)
    {
        buf[i] = rand() % 25;
    }
    return buf;
}

int main (int argc, char** argv)
{
    int number;
    // check for correct number of of arguments
    if (argc != 2) 
    {
        printf("Wrong number of arguments!\n");
        return EXIT_FAILURE;
    }
    // check if argument is integer
    if(!sscanf(argv[1], "%i", &number))
    {
        printf("Could not parse argument!\n");
        return EXIT_FAILURE;
    }

    
    // init mpi
    MPI_Init(&argc, &argv);

    // array length
    int N = number;
    
    // the buffer
    int* buf;
    // first element to know when we are done
    int firstElement;
    // get number of processes
   	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    // get rank of this process
	 	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // check if array ist larger than number of processes
    if(nprocs > N)
    {
      if (rank == 0) {
        printf("Too many processes for the defined array length\n");
      }
      // finalize all processes
      MPI_Finalize();
      return EXIT_FAILURE;
    }
    // array length divided by processes
    int size = N / nprocs;
    // 
    int sizeRemainder = N % nprocs;
    // size of each buffer
    int bufSize;
    // maximum buffer size
    int maxBufSize = size + 1;
    // set right bufSize
    if(rank < sizeRemainder)
    {
        bufSize = maxBufSize;
    }
    else
    {
        bufSize = size;
    }
    // init buffer with right size
    buf = initBuf(bufSize);
    
    // checks that only process with rank 0 prints
    if (rank == 0)
    {
        printf("Vorher: ");
    }
    
    // prints all bufs
    printAllBuffs(buf, bufSize);

    // set first element for finish
    firstElement = buf[0];
    // send firstElement from process with rank 0 to all the others
    MPI_Bcast(&firstElement, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // start circleing
    circle(firstElement, buf, bufSize);
    
    // finalize all processes
    MPI_Finalize();

    // prettier in terminal output
    printf("\n"); 

    return EXIT_SUCCESS;
}