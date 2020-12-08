/****************************************************************************/
/****************************************************************************/
/**                                                                        **/
/**                 TU München - Institut für Informatik                   **/
/**                                                                        **/
/** Copyright: Prof. Dr. Thomas Ludwig                                     **/
/**            Andreas C. Schmidt                                          **/
/**                                                                        **/
/** File:      partdiff.c                                                  **/
/**                                                                        **/
/** Purpose:   Partial differential equation solver for Gauß-Seidel and    **/
/**            Jacobi method.                                              **/
/**                                                                        **/
/****************************************************************************/
/****************************************************************************/

/* ************************************************************************ */
/* Include standard header file.                                            */
/* ************************************************************************ */
#define _POSIX_C_SOURCE 200809L

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <malloc.h>
#include <sys/time.h>
#include <pthread.h>

#include "partdiff.h"

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

struct thread_arguments
{
	double** Matrix_In; // global Matrix_In
	double** Matrix_Out; // global Matrix_out
	double pih; // global value of pih
	double fpisin; // global value of fpisin
	int start; // the row the thread should start at
	int end;// the row the thread should end at
	int N; // global N
	struct options* options; // global options
  int* term_iteration; // pointer to current term_iteration
	double* maxResiduum; // pointer to current maxResiduum
	pthread_mutex_t* maxResiduumMutex; // pointer to maxResiduum Mutex to lock when reading and writing to maxResiduum
};

/* ************************************************************************ */
/* Global variables                                                         */
/* ************************************************************************ */

/* time measurement variables */
struct timeval start_time; /* time when program started                      */
struct timeval comp_time;  /* time when calculation completed                */


/* ************************************************************************ */
/* initVariables: Initializes some global variables                         */
/* ************************************************************************ */
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

/* ************************************************************************ */
/* freeMatrices: frees memory for matrices                                  */
/* ************************************************************************ */
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

/* ************************************************************************ */
/* allocateMemory ()                                                        */
/* allocates memory and quits if there was a memory allocation problem      */
/* ************************************************************************ */
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

/* ************************************************************************ */
/* allocateMatrices: allocates memory for matrices                          */
/* ************************************************************************ */
static
void
allocateMatrices (struct calculation_arguments* arguments)
{
	uint64_t i, j;

	uint64_t const N = arguments->N;

	arguments->M = allocateMemory(arguments->num_matrices * (N + 1) * (N + 1) * sizeof(double));
	arguments->Matrix = allocateMemory(arguments->num_matrices * sizeof(double**));

	for (i = 0; i < arguments->num_matrices; i++)
	{
		arguments->Matrix[i] = allocateMemory((N + 1) * sizeof(double*));

		for (j = 0; j <= N; j++)
		{
			arguments->Matrix[i][j] = arguments->M + (i * (N + 1) * (N + 1)) + (j * (N + 1));
		}
	}
}

/* ************************************************************************ */
/* initMatrices: Initialize matrix/matrices and some global variables       */
/* ************************************************************************ */
static
void
initMatrices (struct calculation_arguments* arguments, struct options const* options)
{
	uint64_t g, i, j; /* local variables for loops */

	uint64_t const N = arguments->N;
	double const h = arguments->h;
	double*** Matrix = arguments->Matrix;

	/* initialize matrix/matrices with zeros */
	for (g = 0; g < arguments->num_matrices; g++)
	{
		for (i = 0; i <= N; i++)
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
			for (i = 0; i <= N; i++)
			{
				Matrix[g][i][0] = 1.0 - (h * i);
				Matrix[g][i][N] = h * i;
				Matrix[g][0][i] = 1.0 - (h * i);
				Matrix[g][N][i] = h * i;
			}

			Matrix[g][N][0] = 0.0;
			Matrix[g][0][N] = 0.0;
		}
	}
}


/* ************************************************************************ */
/* calculate: calculates for specified range of rows per thread             */
/* ************************************************************************ */
static
void*
thread_calculate (void* args)
{
	int i, j;           /* local variable for inner loop */
	double star;        /* four times center value minus 4 neigh.b values */
	double residuum;    /* residuum of current iteration */
	struct thread_arguments *t_args = (struct thread_arguments*) args; // the thread_arguments for this this call
	const struct options *options = t_args->options; // the global options

  double private_max_res = 0; // compute maxResiduum on every Thread and compare after they're done

	/* over all rows */
	for (i = t_args->start; i < t_args->end; i++) // i starts at thread_arguments start and ends at thread_arguments end
	{
		double fpisin_i = 0.0;

		if (options->inf_func == FUNC_FPISIN)
		{
			fpisin_i = t_args->fpisin * sin(t_args->pih * (double)i); // fpisin of thread_arguments
		}

		/* over all columns */
		for (j = 1; j < t_args->N; j++)
		{
			star = 0.25 * (t_args->Matrix_In[i-1][j] + t_args->Matrix_In[i][j-1] 
      + t_args->Matrix_In[i][j+1] + t_args->Matrix_In[i+1][j]); // all Matrices from thread_arguments
			if (options->inf_func == FUNC_FPISIN) // options from thread_arguments 
			{
				star += fpisin_i * sin(t_args->pih * (double)j); // fpisin of thread_arguments
			}
			if (options->termination == TERM_PREC || *(t_args->term_iteration) == 1)
			{
				residuum = t_args->Matrix_In[i][j] - star;
				residuum = (residuum < 0) ? -residuum : residuum;
				private_max_res = (residuum < private_max_res) ? private_max_res : residuum; // update private maxResiduum
			}
			t_args->Matrix_Out[i][j] = star; //Matrix_Out from thread_arguments
		}
	}
  // compare local maxResiduum to global maxResiduum
  pthread_mutex_lock(t_args->maxResiduumMutex); // lock access to maxResiduum
	*(t_args->maxResiduum) = (private_max_res < *(t_args->maxResiduum)) ? 
  *(t_args->maxResiduum) : private_max_res; // update global maxResiduum
	pthread_mutex_unlock(t_args->maxResiduumMutex); // unlock access to maxResiduum
	return NULL;
}

/* ************************************************************************ */
/* calculate: solves the equation                                           */
/* ************************************************************************ */
static
void
calculate (struct calculation_arguments const* arguments, struct calculation_results* results, struct  options * options)
{
	int m1, m2;         /* used as indices for old and new matrices */
	
	double maxResiduum; /* maximum residuum value of a slave in iteration */

	int const N = arguments->N;
	double const h = arguments->h;

	double pih = 0.0;
	double fpisin = 0.0;

	int term_iteration = options->term_iteration;

	pthread_t threads[options->number]; // Array of thread ids 
	pthread_mutex_t maxResiduumMutex; // maxResiduum lock
	int number_of_threads = options->number;
	struct thread_arguments thread_args[number_of_threads]; // Array of thread_arguments to use
	int rows_per_thread = (N-1)/number_of_threads; // nuber of rows each thread should calculate the solution for

	pthread_mutex_init(&maxResiduumMutex, NULL); // init the lock

	/* initialize m1 and m2 depending on algorithm */
	if (options->method == METH_JACOBI)
	{
		m1 = 0;
		m2 = 1;
	}
	else
	{
		m1 = 0;
		m2 = 0;
	}

	if (options->inf_func == FUNC_FPISIN)
	{
		pih = PI * h;
		fpisin = 0.25 * TWO_PI_SQUARE * h * h;
	}

	for (int i = 0; i < number_of_threads; i++)
	{
		thread_args[i].term_iteration = &term_iteration;
		thread_args[i].pih = pih;
		thread_args[i].fpisin = fpisin;
		thread_args[i].N = N;
		thread_args[i].options = options;
		thread_args[i].maxResiduum = &maxResiduum;
		thread_args[i].maxResiduumMutex = &maxResiduumMutex;
	}

	while (term_iteration > 0)
	{
		maxResiduum = 0;
    // loop to asign each thread its rows
		for (int i = 0; i < number_of_threads; i++)
		{
			thread_args[i].start = i * rows_per_thread + 1; 
			if (thread_args[i].end > N-1)
			{
				thread_args[i].end = N-1;
			}
			else
			{
				thread_args[i].end = (i + 1) * rows_per_thread + 1;
			}
			thread_args[i].Matrix_Out = arguments->Matrix[m1];
			thread_args[i].Matrix_In = arguments->Matrix[m2];
			pthread_create(&threads[i], NULL, &thread_calculate, &thread_args[i]);
		}
    //loop to wait for each thread to finish
		for(int k = 0; k < number_of_threads; k++)
		{
			pthread_join(threads[k], NULL);
		}

		results->stat_iteration++;
		results->stat_precision = maxResiduum;

		/* exchange m1 and m2 */
		int j = m1;
		m1 = m2;
		m2 = j;

		/* check for stopping calculation depending on termination method */
		if (options->termination == TERM_PREC)
		{
			if (maxResiduum < options->term_precision)
			{
				term_iteration = 0;
			}
		}
		else if (options->termination == TERM_ITER)
		{
			term_iteration--;
		}
	}
	results->m = m2;
}

/* ************************************************************************ */
/*  displayStatistics: displays some statistics about the calculation       */
/* ************************************************************************ */
static
void
displayStatistics (struct calculation_arguments const* arguments, struct calculation_results const* results, struct options const* options)
{
	int N = arguments->N;
	double time = (comp_time.tv_sec - start_time.tv_sec) + (comp_time.tv_usec - start_time.tv_usec) * 1e-6;

	printf("Berechnungszeit:    %f s \n", time);
	printf("Speicherbedarf:     %f MiB\n", (N + 1) * (N + 1) * sizeof(double) * arguments->num_matrices / 1024.0 / 1024.0);
	printf("Berechnungsmethode: ");

	if (options->method == METH_GAUSS_SEIDEL)
	{
		printf("Gauß-Seidel");
	}
	else if (options->method == METH_JACOBI)
	{
		printf("Jacobi");
	}

	printf("\n");
	printf("Interlines:         %" PRIu64 "\n",options->interlines);
	printf("Stoerfunktion:      ");

	if (options->inf_func == FUNC_F0)
	{
		printf("f(x,y) = 0");
	}
	else if (options->inf_func == FUNC_FPISIN)
	{
		printf("f(x,y) = 2pi^2*sin(pi*x)sin(pi*y)");
	}

	printf("\n");
	printf("Terminierung:       ");

	if (options->termination == TERM_PREC)
	{
		printf("Hinreichende Genaugkeit");
	}
	else if (options->termination == TERM_ITER)
	{
		printf("Anzahl der Iterationen");
	}

	printf("\n");
	printf("Anzahl Iterationen: %" PRIu64 "\n", results->stat_iteration);
	printf("Norm des Fehlers:   %e\n", results->stat_precision);
	printf("\n");
}

/****************************************************************************/
/** Beschreibung der Funktion displayMatrix:                               **/
/**                                                                        **/
/** Die Funktion displayMatrix gibt eine Matrix                            **/
/** in einer "ubersichtlichen Art und Weise auf die Standardausgabe aus.   **/
/**                                                                        **/
/** Die "Ubersichtlichkeit wird erreicht, indem nur ein Teil der Matrix    **/
/** ausgegeben wird. Aus der Matrix werden die Randzeilen/-spalten sowie   **/
/** sieben Zwischenzeilen ausgegeben.                                      **/
/****************************************************************************/
static
void
displayMatrix (struct calculation_arguments* arguments, struct calculation_results* results, struct options* options)
{
	int x, y;

	double** Matrix = arguments->Matrix[results->m];

	int const interlines = options->interlines;

	printf("Matrix:\n");

	for (y = 0; y < 9; y++)
	{
		for (x = 0; x < 9; x++)
		{
			printf ("%7.4f", Matrix[y * (interlines + 1)][x * (interlines + 1)]);
		}

		printf ("\n");
	}

	fflush (stdout);
}

/* ************************************************************************ */
/*  main                                                                    */
/* ************************************************************************ */
int
main (int argc, char** argv)
{
	struct options options;
	struct calculation_arguments arguments;
	struct calculation_results results;

	askParams(&options, argc, argv);

	initVariables(&arguments, &results, &options);

	allocateMatrices(&arguments);
	initMatrices(&arguments, &options);

	gettimeofday(&start_time, NULL);
	calculate(&arguments, &results, &options);
	gettimeofday(&comp_time, NULL);

	displayStatistics(&arguments, &results, &options);
	displayMatrix(&arguments, &results, &options);

	freeMatrices(&arguments);

	return 0;
}
