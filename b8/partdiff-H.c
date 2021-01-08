#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>

#include <mpi.h>


#include "partdiff.h"

// rank of this processs and number of total processes
int rank, nprocs;

struct calculation_arguments
{
	uint64_t  N;              /* number of spaces between lines (lines=N+1)     */
	uint64_t  num_matrices;   /* number of matrices                             */
	double    h;              /* length of a space between two lines            */
	double    ***Matrix;      /* index matrix used for addressing M             */
	double    *M;             /* two matrices with real values                  */
};

struct calculation_results
{
	uint64_t  m;
	uint64_t  stat_iteration; /* number of current iteration                    */
	double    stat_precision; /* actual precision of all slaves in iteration    */
};

static
void*
allocateMemory (size_t size)
{
	void *p;

	if ((p = malloc(size)) == NULL)
	{
		printf("Speicherprobleme! (%" PRIu64 " Bytes angefordert)\n", size);
		exit(1);
	}

	return p;
}
/**
 * rank and size are the MPI rank and size, respectively.
 * from and to denote the global(!) range of lines that this process is responsible for.
 *
 * Example with 9 matrix lines and 4 processes:
 * - rank 0 is responsible for 1-2, rank 1 for 3-4, rank 2 for 5-6 and rank 3 for 7.
 *   Lines 0 and 8 are not included because they are not calculated.
 * - Each process stores two halo lines in its matrix (except for ranks 0 and 3 that only store one).
 * - For instance: Rank 2 has four lines 0-3 but only calculates 1-2 because 0 and 3 are halo lines for other processes. It is responsible for (global) lines 5-6.
 */
static
void
DisplayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options, int rank, int size, int from, int to)
{
  int const elements = 8 * options->interlines + 9;

  int x, y;
  double** Matrix = arguments->Matrix[results->m];
  MPI_Status status;

  /* first line belongs to rank 0 */
  if (rank == 0)
    from--;

  /* last line belongs to rank size - 1 */
  if (rank + 1 == size)
    to++;

  if (rank == 0)
    printf("Matrix:\n");

  for (y = 0; y < 9; y++)
  {
    int line = y * (options->interlines + 1);

    if (rank == 0)
    {
      /* check whether this line belongs to rank 0 */
      if (line < from || line > to)
      {
        /* use the tag to receive the lines in the correct order
         * the line is stored in Matrix[0], because we do not need it anymore */
        MPI_Recv(Matrix[0], elements, MPI_DOUBLE, MPI_ANY_SOURCE, 42 + y, MPI_COMM_WORLD, &status);
      }
    }
    else
    {
      if (line >= from && line <= to)
      {
        /* if the line belongs to this process, send it to rank 0
         * (line - from + 1) is used to calculate the correct local address */
        MPI_Send(Matrix[line - from + 1], elements, MPI_DOUBLE, 0, 42 + y, MPI_COMM_WORLD);
      }
    }

    if (rank == 0)
    {
      for (x = 0; x < 9; x++)
      {
        int col = x * (options->interlines + 1);

        if (line >= from && line <= to)
        {
          /* this line belongs to rank 0 */
          printf("%7.4f", Matrix[line][col]);
        }
        else
        {
          /* this line belongs to another rank and was received above */
          printf("%7.4f", Matrix[0][col]);
        }
      }

      printf("\n");
    }
  }

  fflush(stdout);
}

static
void
alloctMatrices (struct calculation_arguments* arguments,int from, int to) {
  int rows = to - from;
  //int const elements = 8 * options->interlines + 9;
  uint64_t i, j;
  
	uint64_t const N = arguments->N;
  // size of the rows*cols
  int rowssize = (N+1) * rows * sizeof(double);

  arguments->M = allocateMemory(arguments->num_matrices * rowssize);
  arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

  for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory(rowssize * sizeof(double*));
		for (j = 0; j <= N; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * rowssize) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options,int from, int to)
{
  uint64_t rows = (uint64_t)to - (uint64_t)from;

	uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= rows; i++)
		{
			for (j = 0; j <= N; j++)
			{
				Matrix[g][i][j] = 0.0;
			}
		}
	}

	/* initialize borders, depending on function (function 2: nothing to do) */
	if (options->inf_func == FUNC_F0)
	{
		for (g = 0; g < arguments->num_matrices; g++)
		{
      for (i = 0; i <= rows; i++)
			{
			  Matrix[g][i][0] = 1.0 - (h * (i + from));
			  Matrix[g][i][N] = h * (i + from);
			}
      // init first col
      if (rank == 0)
      {
        //set index 0 = 1.0 - (h * i)
        for (i = 0; i <= N; i++)
			  {
          Matrix[g][0][i] = 1.0 - (h * i);
			  }
        Matrix[g][0][N] = 0;
      }
      // init las col
      else if (rank + 1 == nprocs) 
      {
        // set last col to 1 -(h*i)
        for (i = 0; i <= N; i++)
			  {
          Matrix[g][rows][i] = h * i;
			  }
        Matrix[g][N][0] = 0;
      }
		}
	}
}

static
void
initVariables (struct calculation_arguments* arguments, struct calculation_results* results, struct options const* options)
{
	arguments->N = (options->interlines * 8) + 9 - 1;
	arguments->num_matrices = (options->method == METH_JACOBI) ? 2 : 1;
	arguments->h = 1.0 / arguments->N;

	results->m = 0;
	results->stat_iteration = 0;
	results->stat_precision = 0;
}



static
void
freeMatrices (struct calculation_arguments* arguments)
{
	uint64_t i;

	for (i = 0; i < arguments->num_matrices; i++)
	{
		free(arguments->Matrix[i]);
	}

	free(arguments->Matrix);
	free(arguments->M);
}

int main(int argc, char** argv) {

  struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

  askParams(&options, argc, argv);
	initVariables(&arguments, &results, &options);

  uint64_t const N = arguments.N;

  // init mpi
  MPI_Init(&argc, &argv);
  
  // get nr of processes
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  // get rank of this process
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  int rowsPerRank = (N + 1) / nprocs;
  int from = rank * rowsPerRank;
  int to;
  // last rank gets all the rows that arent assigned yet
  if (rank == nprocs - 1)
  {
    to = N;
  }
  else
  {
    to = from + rowsPerRank;
  }

  alloctMatrices(&arguments,from,to);
  initMatrices(&arguments, &options,from,to);

  DisplayMatrix(&arguments,&results, &options, rank, nprocs, from, to);

  MPI_Finalize();

  return EXIT_SUCCESS;
}